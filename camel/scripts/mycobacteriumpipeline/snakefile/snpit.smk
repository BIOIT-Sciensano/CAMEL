from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import variant_calling


rule snpit_report:
    """
    Creates the report for the snpit assay.
    """
    input:
        VCF = variant_calling.get_vcf(config)
    output:
        VAL_HTML = 'snpit/html.iob', # snpit.OUTPUT_REPORT
        INFORMS = 'snpit/informs.io' # snpit.OUTPUT_INFORMS
    params:
        dir_ = 'snpit'
    run:
        from camel.app.tools.snpit.snpit import Snpit
        snpit = Snpit()
        snakemakeutils.add_pickle_inputs(snpit, input)
        step = Step(rule_name=str(rule), tool=snpit, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(snpit, output)

rule snpit_report_empty:
    """
    Creates an empty report when the snpit assay is disabled.
    """
    output:
        HTML = 'snpit/html-empty.iob' # snpit.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Snpit', Path(output.HTML), 2)

rule snpit_export_summary:
    """
    Exports the summary information for the snpit assay.
    """
    input:
        INFORMS = rules.snpit_report.output.INFORMS
    output:
        FILE = 'snpit/summary_out.{ext}' # snpit.OUTPUT_SNPIT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        snpit_informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_out = [(f'snpit_{k}', snpit_informs[k]) for k in ('species', 'lineage', 'sublineage', 'percent_matched')]
        snakemakeutils.export_summary(data_out, Path(output.FILE), str(params.ext), 'snpit')
