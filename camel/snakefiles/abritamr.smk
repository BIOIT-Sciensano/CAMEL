from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly


rule abritamr_run:
    """
    This rule executes AbriTAMR run and gets the results.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TXT_matches ='abritamr/run/txt_matches.io',
        TXT_partials ='abritamr/run/txt_partials.io',
        INFORMS = 'abritamr/run/informs.io'
    params:
        db_path_amrf = config['abritamr']['amrfinderplus']['path'],
        species = config['abritamr']['species']
    run:
        from camel.app.tools.abritamr.abritamrrun import AbriTAMRRun

        abritamr_run_tool = AbriTAMRRun()
        abritamr_run_tool.add_input_files({'DIR_AMRF': [ToolIODirectory(Path(str(params.db_path_amrf)))]})
        snakemakeutils.add_io_input(abritamr_run_tool,'FASTA', Path(input.FASTA))
        step = Step(rule_name=str(rule), tool=abritamr_run_tool, dir_=snakemakeutils.get_rule_dir(output))
        abritamr_run_tool.update_parameters(species=params.species)
        step.run()
        snakemakeutils.dump_io_outputs(abritamr_run_tool, output)

rule abritamr_generate_dummy_mdu_qc_file:
    """
    This rule creates a dummy QC file required by AbriTAMR report in order to make it work as expected.
    Note: The directory specified in the QC file must exactly match the output directory of the AbriTAMR run.
    """
    input:
        TXT_partials = rules.abritamr_run.output.TXT_partials
    output:
        TXT = 'abritamr/dummy_qc/txt.io'
    params:
        species = config['abritamr']['species'],
        dir_ = 'abritamr/dummy_qc'
    run:
        # Extract the working directory of abritAMR
        path_partials = snakemakeutils.load_object(Path(input.TXT_partials))[0].path
        wd_partials = str(path_partials.parent)

        # Create the dummy QC file
        species = params.species if not params.species == 'Salmonella' else 'Salmonella enterica'
        qc_file_txt = Path(params.dir_, 'qc_file.txt')
        with qc_file_txt.open('w') as handle:
            handle.write(f'ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n{wd_partials},{species},{species},PASS\n')
        snakemakeutils.dump_object([ToolIOFile(qc_file_txt)], Path(output.TXT))

rule abritamr_report_run:
    """
    This rule will run AbriTAMR report to generate an antibiogram in case the species is Salmonella,
    else it will create a summary of the TXT_MATCHES and TXT_PARTIALS that will not actually be used in the
    abritamr_report rule.
    """
    input:
        TXT_mdu_qc = rules.abritamr_generate_dummy_mdu_qc_file.output.TXT,
        TXT_matches = rules.abritamr_run.output.TXT_matches,
        TXT_partials = rules.abritamr_run.output.TXT_partials,
        INFORMS_abritamr_run = rules.abritamr_run.output.INFORMS
    output:
        REPORT_abritamr = 'abritamr/report/report.io',
        INFORMS = 'abritamr/report/informs.io'
    run:
        from camel.app.tools.abritamr.abritamrreport import AbriTAMRReport
        abritamr_report_tool = AbriTAMRReport()
        snakemakeutils.add_io_inputs(abritamr_report_tool, input)
        step = Step(rule_name=str(rule), tool=abritamr_report_tool, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(abritamr_report_tool, output)

rule abritamr_create_summary:
    """
    Creates the tabular summary output for the AbritAMR assay.
    """
    input:
        REPORT_abritamr = rules.abritamr_report_run.output.REPORT_abritamr,
        INFORMS_abritamr_run = rules.abritamr_run.output.INFORMS,
        INFORMS_abritamr_report = rules.abritamr_report_run.output.INFORMS
    output:
        TSV = 'abritamr/summary/summary_out.{ext}' # abritamr.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import pandas as pd

        df_abritamr = pd.read_excel(
            snakemakeutils.load_object(Path(input.REPORT_abritamr))[0].path, engine='openpyxl', dtype=str)
        df_abritamr.fillna('-', inplace=True)  # replace all missing values by dashes
        data_summary = []
        with Path(output.TSV).open('w') as handle:
            for column in df_abritamr.columns[2:]:
                key = f'abritamr_{column.replace(" - ","_")}'
                value = df_abritamr.iloc[0][column]
                data_summary.append((key, value))

            informs_abritamr_run = snakemakeutils.load_object(Path(input.INFORMS_abritamr_run))
            informs_abritamr_report = snakemakeutils.load_object(Path(input.INFORMS_abritamr_report))
            data_summary.append(('abritamr_tool_version', informs_abritamr_report['_version']))
            data_summary.append(('abritamr_db_version', informs_abritamr_run['last_update_date']))
        snakemakeutils.export_summary(data_summary, Path(output.TSV), str(params.ext), 'abritAMR')

rule abritamr_report:
    """
    This rule creates a simple HTML output report for the AbriTAMR tool.
    """
    input:
        REPORT_abritamr = rules.abritamr_report_run.output.REPORT_abritamr,
        INFORMS_abritamr_run = rules.abritamr_run.output.INFORMS,
        TXT_matches = rules.abritamr_run.output.TXT_matches,
        TXT_partials = rules.abritamr_run.output.TXT_partials
    output:
        VAL_HTML = 'abritamr/output_report/html.iob' # abritamr.OUTPUT_REPORT
    run:
        from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter

        reporter = AbriTAMRReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule abritamr_report_empty:
    """
    Creates an empty HTML report for the AbriTAMR analysis.
    """
    output:
        VAL_HTML = 'abritamr/output_report/html-empty.iob' # abritamr.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.abritamr.abritamrreporter import AbriTAMRReporter
        snakepipelineutils.create_empty_report_section(AbriTAMRReporter.TITLE, Path(output.VAL_HTML))
