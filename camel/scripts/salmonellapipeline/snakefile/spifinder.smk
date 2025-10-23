from pathlib import Path

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.snakefiles import assembly
from camel.scripts.salmonellapipeline.snakefile import spifinder


rule spifinder_fastq_run:
    """
    Runs SPIFinder on FASTQ input. 
    """
    input:
        IO = 'fq_dict.io'
    output:
        JSON = 'spifinder/spifinder_fastq/json.io', # spifinder.OUTPUT_FASTQ_JSON
        INFORMS = 'spifinder/spifinder_fastq/informs.io' # spifinder.OUTPUT_FASTQ_INFORMS
    params:
        dir_ = 'spifinder/spifinder_fastq',
        db_path = config['spifinder']['path']
    run:
        from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
        spifinder_tool = SPIFinder()
        spifinder_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db_path))]})
        if config['input_type'] == 'illumina':
            spifinder_tool.add_input_files(snakepipelineutils.extract_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        elif config['input_type'] == 'ont':
            spifinder_tool.add_input_files(snakepipelineutils.extract_fq_input(
                Path(input.IO), key_se='FASTQ', read_type='SE'))
        step = Step(rule_name=str(rule), tool=spifinder_tool, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(spifinder_tool, output)

rule spifinder_fasta_run:
    """
    Runs SPIFinder on FASTA input.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        JSON = 'spifinder/fasta/json.io', # spifinder.OUTPUT_FASTA_JSON
        INFORMS = 'spifinder/fasta/informs.io' # spifinder.OUTPUT_FASTA_INFORMS
    params:
        dir_ = 'spifinder/fasta',
        db_path = config['spifinder']['path']
    run:
        from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
        spifinder_tool = SPIFinder()
        spifinder_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db_path))]})
        snakemakeutils.add_pickle_input(spifinder_tool, 'FASTA', Path(input.FASTA))
        step = Step(rule_name=str(rule), tool=spifinder_tool, dir_=Path(params.dir_))
        step.run()
        spifinder_tool.informs['_tag'] = 'FASTA'
        snakemakeutils.dump_tool_outputs(spifinder_tool, output)

rule spifinder_create_summary:
    """
    This rule creates a summary output for the hits of SPIFinder in fastq and fasta mode.
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if config['input_type'] in ('ont', 'illumina') else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS
    output:
        FILE = 'spifinder/summary/summary_out.{ext}' # spifinder.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs_fasta = snakemakeutils.load_object(Path(input.INFORMS_spifinder_fasta))
        data_summary = [
            ('spifinder_tool_version', informs_fasta['_name_full']),
            ('spifinder_db_version', informs_fasta['last_update_date']),
        ]

        # Parse FASTA
        path_in = snakemakeutils.load_object(Path(input.JSON_FASTA))[0].path
        data_fasta = spifinder.parse_json(path_in,'fasta')

        # Parse FASTQ
        if len(input.JSON_FASTQ) > 0:
            path_in = snakemakeutils.load_object(Path(input.JSON_FASTQ))[0].path
            data_fastq = spifinder.parse_json(path_in,'fastq')
        else:
            data_fastq = None

        # Construct output
        if params.ext == 'tsv':
            data_summary.append(('spifinder_fasta', str(data_fasta.values.tolist())))
            if data_fastq is not None:
                data_summary.append(('spifinder_fastq', str(data_fastq.values.tolist())))
        elif params.ext == 'json':
            data_summary.append(('spifinder_fasta', data_fasta.to_dict('records')))
            if data_fastq is not None:
                data_summary.append(('spifinder_fastq', data_fastq.to_dict('records')))
        else:
            raise ValueError(f'Invalid ext: {params.ext}')
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'spifinder')

rule spifinder_report:
    """
    This rule creates a simple output report, combining both SPIFinder tables in one report.
    """
    input:
        JSON_FASTQ = rules.spifinder_fastq_run.output.JSON if config['input_type'] in ('ont', 'illumina') else [],
        JSON_FASTA = rules.spifinder_fasta_run.output.JSON,
        INFORMS_spifinder_fastq = rules.spifinder_fastq_run.output.INFORMS if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_spifinder_fasta = rules.spifinder_fasta_run.output.INFORMS,
        CSV_metadata = config['spifinder']['metadata']
    output:
        VAL_HTML = 'spifinder/report/html.iob' # spifinder.OUTPUT_REPORT
    params:
        dir_ = 'spifinder/report'
    run:
        import pandas as pd
        from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter

        # Export documentation
        file = pd.read_csv(input.CSV_metadata, delimiter=';')
        path_doc = Path(params.dir_, 'documentation.tsv')
        file.to_csv(path_doc, sep='\t')

        # Create the report
        reporter = SPIFinderReporter()
        snakemakeutils.add_pickle_inputs(reporter, input, excluded_keys=['JSON_FASTQ', 'INFORMS_spifinder_fastq', 'CSV_metadata'])
        reporter.add_input_files({'TSV_documentation': [ToolIOFile(path_doc)]})
        if input.JSON_FASTQ:
            snakemakeutils.add_pickle_input(reporter, 'JSON_FASTQ', Path(input.JSON_FASTQ))
        if input.INFORMS_spifinder_fastq:
            reporter.add_input_informs({'spifinder_fastq': snakemakeutils.load_object(Path(input.INFORMS_spifinder_fastq))})
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_output(reporter, 'VAL_HTML', Path(output.VAL_HTML))
        # snakemakeutils.dump_object([ToolIOFile(path_doc)], Path(output.TSV_doc))

rule spifinder_report_empty:
    """
    Creates an empty HTML report for the SPIFinder analysis.
    """
    output:
        VAL_HTML = 'spifinder/report/html-empty.iob' # # spifinder.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section(SPIFinderReporter.TITLE, Path(output.VAL_HTML))
