from pathlib import Path

from camel.resources.snakefile import trimming_illumina, gene_detection, trimming, variant_calling, variant_filtering, \
    contamination_check_kraken, quality_checks, sequence_typing, downsampling, amrfinder, mobsuite, quast, confindr, \
    core, trimming_ont, assembly, human_read_scrubbing, resfinder4, bacmet, read_simulation
from camel.scripts.klebsiellapipeline.snakefile import kleborate

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE
include: human_read_scrubbing.SNAKEFILE
include: read_simulation.SNAKEFILE
include: downsampling.SNAKEFILE
include: trimming_illumina.SNAKEFILE
include: trimming_ont.SNAKEFILE
include: assembly.SNAKEFILE
include: quast.SNAKEFILE
include: contamination_check_kraken.SNAKEFILE
include: confindr.SNAKEFILE
include: quality_checks.SNAKEFILE
include: variant_calling.SNAKEFILE
include: variant_filtering.SNAKEFILE
include: gene_detection.SNAKEFILE
include: sequence_typing.SNAKEFILE
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: kleborate.SNAKEFILE
include: mobsuite.SNAKEFILE
include: bacmet.SNAKEFILE


#########
# Rules #
#########
rule all:
    """
    Rule to generate the required output files.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']

rule report_command_section:
    """
    Creates a report section with the commands used in the pipeline. 
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = quast.OUTPUT_INFORMS,
        INFORMS_busco = quast.OUTPUT_INFORMS_BUSCO,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        # INFORMS_mapping = quality_checks.get_mapping_rate_informs(config),
        # INFORMS_depth = quality_checks.get_depth_informs(config),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        INFORMS_amrfnder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS  if 'mob_suite' in config['analyses'] else [],
        INFORMS_kleborate = kleborate.OUTPUT_INFORMS if 'kleborate' in config['analyses'] else [],
        INFORMS_vfdb_core = str(gene_detection.OUTPUT_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        IFNORMS_bacmet = bacmet.OUTPUT_INFORMS if 'bacmet' in config['analyses'] else [],
        IFNORMS_prodigal = bacmet.OUTPUT_PRODIGAL_INFORMS if 'bacmet' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = quast.OUTPUT_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = quality_checks.OUTPUT_REPORT.format(input_type=config['input_type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # AMR detection
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_resfinder4 = resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY,
        # Virulence gene detection
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = 'mob_suite/genomic_context/html.iob',
        # Kleborate
        report_kleborate = kleborate.OUTPUT_REPORT if 'kleborate' in config['analyses'] else kleborate.OUTPUT_REPORT_EMPTY,
        # BacMet
        report_prodigal = bacmet.OUTPUT_PRODIGAL_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_PRODIGAL_REPORT_EMPTY,
        report_bacmet = bacmet.OUTPUT_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_REPORT_EMPTY,
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_scgmlst = sequence_typing.get_sequence_typing_report('scgmlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir= config['output_dir'],
        pipeline_info = config['pipeline'],
        input_dict = config['input'],
        input_type = config['input_type'],
        detection_method = config['detection_method']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils

        # Add the header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            detection_method=params.detection_method
        ))

        # Set up the report content structure
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure,params.input_type,input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('AMR detection', 'amr', [Path(input.report_amrfinder), Path(input.report_resfinder4)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('Kleborate', 'kleborate', [Path(input.report_kleborate)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Biocide and metal resistance', 'bacmet', [Path(input.report_prodigal), Path(input.report_bacmet)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_scgmlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    Combines the summary information of several steps into a single TSV file.
    """
    input:
        core.OUTPUT_SUMMARY_INIT,
        lambda wildcards: human_read_scrubbing.get_summaries(config, wildcards.ext),
        lambda wildcards: downsampling.get_summaries(config, wildcards.ext),
        trimming.get_summaries(config),
        quast.OUTPUT_SUMMARY,
        quality_checks.OUTPUT_SUMMARY,
        lambda wildcards: contamination_check_kraken.get_summaries(config, wildcards.ext),
        confindr.get_summary(config),
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses'] else [],
        # AMR detection
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses'] else [],
        resfinder4.OUTPUT_SUMMARY if 'resfinder4' in config['analyses'] else [],
        # Virulence detection
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='vfdb_core', ext=wildcards.ext) if 'vfdb_core' in config['analyses'] else [],
        # Kleborate
        kleborate.OUTPUT_SUMMARY if 'kleborate' in config['analyses'] else [],
        # Plasmid characterization
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='plasmidfinder', ext=wildcards.ext) if 'plasmidfinder' in config['analyses'] else [],
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # BacMet
        bacmet.OUTPUT_SUMMARY if 'bacmet' in config['analyses'] else [],
        # Sequence typing
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='scgmlst', ext=wildcards.ext) if 'scgmlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        SnakePipelineUtils.combine_summary_data(input, Path(output.FILE), str(params.ext))

rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = 'amrfinder/tool/tsv.io' if 'amrfinder' in config['analyses'] else [],
        # Virulence
        TSV_gd_vfdb = 'gene_detection/vfdb_core/metadata/tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = str(gene_detection.OUTPUT_DB_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        # BacMet
        TSV_bacmet = 'bacmet/hit_filtering/tsv.io' if 'bacmet' in config['analyses'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
