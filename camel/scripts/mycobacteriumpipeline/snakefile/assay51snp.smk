from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import variant_calling, variant_filtering


rule assay_51snp_init_db:
    """
    Converts the 51 SNP databases to CAMEL IO pickles. 
    """
    output:
        BED = '51snp/db/bed.io',
        TSV = '51snp/db/tsv.io'
    params:
        bed = config['51snp']['bed'],
        tsv = config['51snp']['profiles']
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        snakemakeutils.dump_object([ToolIOFile(Path(params.bed))], Path(output.BED))
        snakemakeutils.dump_object([ToolIOFile(Path(params.tsv))], Path(output.TSV))

rule assay_51snp_filter_vcf:
    """
    Extracts the SNPs at the 51 SNP positions.
    """
    input:
        VCF_GZ = variant_filtering.OUTPUT_VCF,
        BED = rules.assay_51snp_init_db.output.BED
    output:
        VCF = '51snp/filter_vcf/vcf.io'
    params:
        dir_ = '51snp/filter_vcf'
    run:
        from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
        bcf_filter = BcftoolsFilter()
        snakemakeutils.add_pickle_inputs(bcf_filter, input)
        bcf_filter.update_parameters(output_filename='51_snps.vcf')
        step = Step(rule_name=str(rule), tool=bcf_filter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bcf_filter, output)

rule assay_51snp_detect_info:
    """
    Detects the metadata for the 51SNP assay.
    """
    input:
        BED = rules.assay_51snp_init_db.output.BED,
        TSV = rules.assay_51snp_init_db.output.TSV,
        VCF = variant_calling.get_vcf(config),
        VCF_filt = variant_filtering.OUTPUT_VCF
    output:
        INFORMS = '51snp/detect/informs.iob'
    params:
        dir_ = '51snp/detect'
    run:
        from camel.app.tools.pipelines.mycobacterium.assay51snpdetector import Assay51SnpDetector
        detector = Assay51SnpDetector()
        snakemakeutils.add_pickle_inputs(detector, input)
        step = Step(rule_name=str(rule), tool=detector, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(detector, output)

rule assay_51snp_report:
    """
    Creates the report for the 51 SNP assay.
    """
    input:
        INFORMS_detection = rules.assay_51snp_detect_info.output.INFORMS
    output:
        VAL_HTML = '51snp/report/html.iob' # assay51snp.OUTPUT_REPORT
    params:
        dir_ = '51snp/report',
        sample_name = config['sample_name']
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.assay51snpreporter import Assay51SnpReporter
        spr = Assay51SnpReporter()
        snakemakeutils.add_pickle_inputs(spr, input)
        spr.add_input_files({'VAL_Sample': [ToolIOValue(params.sample_name)]})
        step = Step(rule_name=str(rule), tool=spr, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(spr, output)

rule assay_51snp_report_empty:
    """
    Creates an empty report for the 51 SNP species identification.
    """
    output:
        VAL_HTML= '51snp/report/html-empty.iob' # assay51snp.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.pipelines.mycobacterium.assay51snpreporter import Assay51SnpReporter
        snakepipelineutils.create_empty_report_section(Assay51SnpReporter.TITLE, Path(output.VAL_HTML), 3)

rule assay_51snp_dump_summary_info:
    """
    Dumps the summary information from the 51SNP workflow.
    """
    input:
        INFORMS = rules.assay_51snp_detect_info.output.INFORMS
    output:
        FILE = '51snp/summary_out.{ext}' # assay51snp.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [
            ('51SNP-positive_control', informs['mtbc_pos_control']),
            ('51SNP-gyrB_group', informs['gyrB_group']),
            ('51SNP-genetic_group', informs['genetic_group']),
            ('51SNP-scg', informs['scg_profile'].scg),
            ('51SNP-st', informs['scg_profile'].st),
            ('51SNP-matching_snps', informs['scg_nb_snps_matched'])
        ]
        for i in range(1, 52):
            key = 'SNP{:02d}'.format(i)
            data_summary.append((f'51SNP-{key}', informs[key]))
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), '51_snp')
