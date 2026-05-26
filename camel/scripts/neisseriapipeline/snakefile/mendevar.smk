from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import sequence_typing


rule mendevar_run:
    """
    Runs the MenDeVAR stand-alone tool.
    """
    input:
        TSV = sequence_typing.OUTPUT_TSV.format(scheme='bast', locus_type='peptide')
    output:
        TSV = 'mendevar/tsv.io',
        INFORMS = 'mendevar/informs.io'
    params:
        db = config.get('mendevar', {}).get('db')
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.neisseria.mendevar import MenDeVAR

        mendevar_ = MenDeVAR()
        step = Step(rule_name=str(rule), tool=mendevar_, dir_=snakemakeutils.get_rule_dir(output))
        snakemakeutils.add_io_inputs(mendevar_, input)

        # Add database
        if params.db is None:
            raise RuntimeError('MenDeVAR database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'MenDeVAR database not found: {path_db}')
        mendevar_.add_input_files({'DB': [ToolIOFile(path_db)]})

        # Run tool
        step.run()
        snakemakeutils.dump_io_outputs(mendevar_, output)

rule mendevar_report:
    """
    Creates an output report for the MenDeVAR analysis.
    """
    input:
        TSV = rules.mendevar_run.output.TSV,
        INFORMS_mendevar = rules.mendevar_run.output.INFORMS
    output:
        HTML = 'mendevar/report/html.iob' # mendevar.OUTPUT_REPORT
    params:
        dir_ = 'mendevar/report'
    run:
        from camel.app.tools.pipelines.neisseria.mendevarreporter import MenDeVARReporter
        reporter = MenDeVARReporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        snakemakeutils.add_io_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule mendevar_report_empty:
    """
    Creates an empty MenDeVAR report when the analysis is disabled.
    """
    output:
        HTML = 'mendevar/report/html-empty.iob' # mendevar.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('MenDeVAR', Path(output.HTML), 3)

rule mendevar_create_summary:
    """
    Creates the tabular summary output for the MenDeVAR assay.
    """
    input:
        INFORMS_mendevar = rules.mendevar_run.output.INFORMS
    output:
        TSV = 'mendevar/summary/summary_mendevar.{ext}' # mendevar.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        # Collect informs
        informs = snakemakeutils.load_object(Path(input.INFORMS_mendevar))

        # Create output files
        rows_out = [
            ('mendevar_bexsero_status', informs['bexsero_status']),
            ('mendevar_trumenba_status', informs['trumenba_status'])
        ]
        snakemakeutils.export_summary(rows_out, Path(output.TSV), str(params.ext), 'mendevar')
