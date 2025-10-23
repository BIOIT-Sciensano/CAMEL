from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.snakefiles import assembly
from camel.app.core.snakemake import snakemakeutils

rule shigeifinder_run:
    """
    Runs the ShigEiFinder assay.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TSV = 'shigeifinder/tool/tsv.io',
        INFORMS = 'shigeifinder/tool/informs.io' # shigeifinder.OUTPUT_INFORMS
    params:
        dir_ = 'shigeifinder/tool'
    run:
        from camel.app.tools.pipelines.shigella.shigeifinder import ShigEiFinder

        shigeifinder_ = ShigEiFinder()
        step = Step(rule_name=str(rule), tool=shigeifinder_, dir_=Path(str(params.dir_)))
        snakemakeutils.add_pickle_inputs(shigeifinder_,input)

        # Run tool
        step.run()
        snakemakeutils.dump_tool_outputs(shigeifinder_, output)

rule shigeifinder_report:
    """
    Creates an output report for the ShigEiFinder analysis.
    """
    input:
        TSV = rules.shigeifinder_run.output.TSV,
        INFORMS_shigeifinder = rules.shigeifinder_run.output.INFORMS
    output:
        HTML = 'shigeifinder/report/html.iob' # shigeifinder.OUTPUT_REPORT
    params:
        dir_ = 'shigeifinder/report'
    run:
        from camel.app.tools.pipelines.shigella.shigeifinderreporter import ShigEiFinderReporter
        reporter = ShigEiFinderReporter()
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule shigeifinder_report_empty:
    """
    Creates an empty ShigEiFinder report when the analysis is disabled.
    """
    output:
        HTML = 'shigeifinder/report/html-empty.iob' # shigeifinder.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('ShigEiFinder', Path(output.HTML), 3)

rule shigeifinder_create_summary:
    """
    Creates the tabular summary output for the ShigEiFinder assay.
    """
    input:
        INFORMS_shigeifinder = rules.shigeifinder_run.output.INFORMS
    output:
        FILE = 'shigeifinder/summary_shigeifinder.{ext}' # shigeifinder.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS_shigeifinder))
        keys = ['species', 'serotype', 'O_antigen', 'H_antigen', 'cluster']
        data_summary = [(f'shigeifinder_{key}', informs[key]) for key in keys]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'shigeifinder')
