from pathlib import Path

from camel.snakefiles import trimming, trimming_illumina, \
    quality_checks, contamination_check_kraken, variant_calling, variant_filtering, gene_detection, sequence_typing, \
    downsampling, confindr, quast, core, assembly, amrfinder, resfinder4, mobsuite, mykrobe, human_read_scrubbing, \
    read_simulation, trimming_ont
from camel.scripts.shigellapipeline.snakefile import shigeifinder, shigatyper

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE
include: human_read_scrubbing.SNAKEFILE
include: downsampling.SNAKEFILE
include: read_simulation.SNAKEFILE
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
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: mobsuite.SNAKEFILE
include: sequence_typing.SNAKEFILE
include: shigeifinder.SNAKEFILE
include: shigatyper.SNAKEFILE
include: mykrobe.SNAKEFILE


#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
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
        INFORMS_simulation = read_simulation.OUTPUT_INFORMS if config['input_type'] == 'fasta' else [],
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = quast.OUTPUT_INFORMS,
        INFORMS_busco = quast.OUTPUT_INFORMS_BUSCO,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type']),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        INFORMS_shigeifinder = shigeifinder.OUTPUT_INFORMS if 'shigeifinder' in config['analyses'] else [],
        INFORMS_shigatyper = shigatyper.OUTPUT_INFORMS if 'shigatyper' in config['analyses'] else[],
        INFORMS_mykrobe = mykrobe.OUTPUT_INFORMS if 'mykrobe' in config['analyses'] else[],
        INFORMS_amrfinder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_virulence = str(gene_detection.OUTPUT_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_virulence_shiga = str(gene_detection.OUTPUT_INFORMS).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_serotype_o = str(gene_detection.OUTPUT_INFORMS).format(db='serotype_o') if 'serotype_o' in config['analyses'] else [],
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_mlst = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_rmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input,Path(output.HTML), Path(params.dir_))

rule combine_reports:
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
        report_adv_qc = quality_checks.OUTPUT_REPORT,
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # # Gene detection
        report_amrfinder = (amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY),
        report_resfinder4 = (resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY),
        report_virulence = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_virulence_shiga = gene_detection.get_gene_detection_report('virulencefinder_shiga', config, 'virulencefinder'),
        report_serotype_o_type= gene_detection.get_gene_detection_report('serotype_o',config,'serotype_o'),
        # Plasmid characterization
        report_mob_suite = (mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_REPORT_EMPTY),
        report_genomic_context = (mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY),
        # # Shigella serotyping
        report_shigeifinder = (shigeifinder.OUTPUT_REPORT if 'shigeifinder' in config['analyses'] else shigeifinder.OUTPUT_REPORT_EMPTY),
        report_shigatyper = (shigatyper.OUTPUT_REPORT if 'shigatyper' in config['analyses'] else shigatyper.OUTPUT_REPORT_EMPTY),
        report_mykrobe = (mykrobe.OUTPUT_REPORT if 'mykrobe' in config['analyses'] else mykrobe.OUTPUT_REPORT_EMPTY),
        # # Sequence typing
        report_mlst_warwick = sequence_typing.get_sequence_typing_report('mlst_warwick', config),
        report_mlst_pasteur = sequence_typing.get_sequence_typing_report('mlst_pasteur', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name=config['sample_name'],
        output_dir=config['output_dir'],
        pipeline_info=config['pipeline'],
        input_dict=config['input'],
        input_type=config['input_type'],
        detection_method=config['gene_detection']['options']['method'],
        citation_keys=config['citations']
    run:
        import datetime
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.scriptutils.reportpipeline import ReportPipeline

        # Add the header section
        report = snakepipelineutils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(snakepipelineutils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            detection_method=params.detection_method,
            key_citation=params.citation_keys['main']
        ))

        # Add output sections
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure,params.input_type,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('<i>Shigella</i> serotyping', 'shigella_typing', [Path(x) for x in (
                input.report_shigeifinder, input.report_shigatyper, input.report_mykrobe)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_amrfinder, input.report_resfinder4)]),
            ('Virulence characterization', 'viru', [Path(x) for x in (
                input.report_virulence, input.report_virulence_shiga)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst_warwick, input.report_mlst_pasteur, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        snakepipelineutils.add_report_content(report, report_structure)

rule combine_summary_files:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        core.OUTPUT_SUMMARY_INIT,
        lambda wildcards: human_read_scrubbing.get_summaries(config, wildcards.ext),
        lambda wildcards: downsampling.get_summaries(config, wildcards.ext),
        trimming.get_summaries(config),
        quast.OUTPUT_SUMMARY,
        lambda wildcards: contamination_check_kraken.get_summaries(config, wildcards.ext),
        confindr.get_summary(config),
        quality_checks.OUTPUT_SUMMARY,
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        # Shigella typing
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='serotype_o', ext=wildcards.ext) if 'serotype_o' in config['analyses'] else [],
        # Shigella typing
        shigeifinder.OUTPUT_SUMMARY if 'shigeifinder' in config['analyses'] else [],
        shigatyper.OUTPUT_SUMMARY if 'shigatyper' in config['analyses'] else [],
        mykrobe.OUTPUT_SUMMARY if 'mykrobe' in config['analyses'] else [],
         # Gene detection
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses'] else [],
        resfinder4.OUTPUT_SUMMARY if 'resfinder4' in config['analyses'] else [],
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='virulencefinder', ext=wildcards.ext) if 'virulencefinder' in config['analyses'] else [],
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='virulencefinder_shiga', ext=wildcards.ext) if 'virulencefinder' in config['analyses'] else [],
        # Plasmid characterization
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # Sequence typing
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst_pasteur', ext=wildcards.ext) if 'mlst_pasteur' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst_warwick', ext=wildcards.ext) if 'mlst_warwick' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.combine_summary_data(input, Path(output.FILE), str(params.ext))

rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = amrfinder.OUTPUT_TSV if 'amrfinder' in config['analyses'] else [],
        # Virulence
        TSV_gd_coli = 'gene_detection/virulencefinder/metadata/tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_coli = str(gene_detection.OUTPUT_DB_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        TSV_gd_shiga = 'gene_detection/virulencefinder_shiga/metadata/tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_shiga = str(gene_detection.OUTPUT_DB_INFORMS).format(db='virulencefinder_shiga') if 'virulencefinder' in config['analyses'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
