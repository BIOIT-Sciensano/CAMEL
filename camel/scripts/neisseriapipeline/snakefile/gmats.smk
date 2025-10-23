from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import sequence_typing


rule gmats_run:
    """
    Runs the gMATS assay.
    """
    input:
        TSV = sequence_typing.OUTPUT_TSV.format(scheme='bast', locus_type='peptide')
    output:
        TSV = 'gmats/tool/tsv.io',
        INFORMS = 'gmats/tool/informs.io'
    params:
        dir_ = 'gmats/tool',
        db = config.get('gmats', {}).get('db')
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.neisseria.gmats import GMats

        gmats_ = GMats()
        step = Step(rule_name=str(rule), tool=gmats_, dir_=Path(str(params.dir_)))
        snakemakeutils.add_pickle_inputs(gmats_,input)

        # Add database
        if params.db is None:
            raise RuntimeError('gMATS database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'gMATS database not found: {path_db}')
        gmats_.add_input_files({'DB': [ToolIOFile(path_db)]})

        # Run tool
        step.run()
        snakemakeutils.dump_tool_outputs(gmats_, output)

rule gmats_report:
    """
    Creates an output report for the gMATS analysis.
    """
    input:
        TSV = rules.gmats_run.output.TSV,
        INFORMS_gmats = rules.gmats_run.output.INFORMS
    output:
        HTML = 'gmats/report/html.iob' # gmats.OUTPUT_REPORT
    params:
        dir_ = 'gmats/report'
    run:
        from camel.app.tools.pipelines.neisseria.gmatsreporter import GMatsReporter
        reporter = GMatsReporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        snakemakeutils.add_pickle_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule gmats_report_empty:
    """
    Creates an empty gMATS report when the analysis is disabled.
    """
    output:
        HTML = 'gmats/report/html-empty.iob' # gmats.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('gMATS', Path(output.HTML), 3)

rule gmats_create_summary:
    """
    Creates the tabular summary output for the gMATS assay.
    """
    input:
        INFORMS_gmats = rules.gmats_run.output.INFORMS
    output:
        FILE = 'gmats/summary/summary_gmats.{ext}' # gmats.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS_gmats))
        data_summary = [('gmats_status', informs['gMATS_status'])]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'gmats')
