from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import variant_calling, variant_filtering


rule snp_lineage_detection:
    """
    Detects the SNP lineage based on the variants detected in the sample.
    """
    input:
        VCF = variant_calling.get_vcf(config),
        VCF_filt = variant_filtering.OUTPUT_VCF
    output:
        INFORMS = 'snp_lineage/tool/informs.iob' # snplineage.OUTPUT_INFORMS
    params:
        dir_ = 'snp_lineage/tool',
        lineage_snp_positions = config['snp_lineage']['bed']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.mycobacterium.snplineagedetector import SNPLineageDetector
        detector = SNPLineageDetector()
        snakemakeutils.add_pickle_inputs(detector, input)
        detector.add_input_files({'BED': [ToolIOFile(Path(params.lineage_snp_positions))]})
        step = Step(rule_name=str(rule), tool=detector, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(detector, output)

rule snp_lineage_report:
    """
    Creates the HTML output report for the SNP lineage report.
    """
    input:
        INFORMS_detection = rules.snp_lineage_detection.output.INFORMS
    output:
        VAL_HTML = 'snp_lineage/report/html.iob' # snplineage.OUTPUT_REPORT
    params:
        dir_ = 'snp_lineage/report',
        lineage_snp_positions = config['snp_lineage']['bed']
    run:
        from camel.app.tools.pipelines.mycobacterium.snplineagereporter import SNPLineageReporter
        lineage_snp = SNPLineageReporter()
        snakemakeutils.add_pickle_inputs(lineage_snp, input)
        step = Step(rule_name=str(rule), tool=lineage_snp, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(lineage_snp, output)

rule snp_lineage_report_empty:
    """
    Creates an empty report when the assay is disabled.
    """
    output:
        VAL_HTML = 'snp_lineage/report/html-empty.io' # snplineage.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.tools.pipelines.mycobacterium.snplineagereporter import SNPLineageReporter
        from camel.app.snakemake.snakepipelineutils import  SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section(SNPLineageReporter.TITLE, Path(output.VAL_HTML), 2)

rule snp_lineage_dump_summary_info:
    """
    Dumps the summary information from the SNP lineage assay.
    """
    input:
        INFORMS = rules.snp_lineage_detection.output.INFORMS
    output:
        FILE = 'snp_lineage/summary/summary_out.{ext}' # snplineage.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        all_lineages = ', '.join([
            d['lineage'].id_ for _, d in informs['detected_lineage_by_level'].items() if d is not None])
        data_summary = [('snp_lineages', all_lineages)]
        if params.ext == 'json' and informs.get('detected_lineage_by_level'):
            for key, value in informs['detected_lineage_by_level'].items():
                if value is not None and value.get('lineage'):
                    informs['detected_lineage_by_level'][key]['lineage'] = value['lineage'].__dict__
            data_summary.append(('detected_lineage_by_level', informs['detected_lineage_by_level']))
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'snplineage')
