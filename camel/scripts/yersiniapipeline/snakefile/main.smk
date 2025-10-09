from pathlib import Path

from camel.resources.snakefile import trimming, trimming_illumina, quality_checks, \
    contamination_check_kraken, gene_detection, sequence_typing, downsampling, confindr, quast, core, trimming_ont, \
    mobsuite, resfinder4, amrfinder, assembly, human_read_scrubbing, read_simulation
from camel.scripts.yersiniapipeline.snakefile import species_determination

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
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: gene_detection.SNAKEFILE
include: mobsuite.SNAKEFILE
include: sequence_typing.SNAKEFILE
include: species_determination.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    This rule ensures that the required output files are generated.
    """
    input:
        config['output_report'],
        config['output_tabular']

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
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
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input_type']),
        INFORMS_amrfinder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_vfdb_core = str(gene_detection.OUTPUT_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_mlst = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_mlst_mcnally = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst_mcnally') if 'mlst_mcnally' in config['analyses'] else [],
        INFORMS_amr = sequence_typing.OUTPUT_INFORMS.format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        INFORMS_ampc = str(gene_detection.OUTPUT_INFORMS).format(db='ampc') if 'ampc' in config['analyses'] else[],
        INFORMS_cgmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_cgmlst_ye = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst_ye') if 'cgmlst_ye' in config['analyses'] else [],
        INFORMS_cgmlst_yp = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst_yp') if 'cgmlst_yp' in config['analyses'] else [],
        INFORMS_cgmlst_enterobase = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst_enterobase') if 'cgmlst_enterobase' in config['analyses'] else [],
        INFORMS_rmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='rmlst') if 'rmlst' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), Path(params.dir_))

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
        report_adv_qc = str(quality_checks.OUTPUT_REPORT).format(
            input_type=config['input_type']),
        # AMR detection
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_resfinder4 = resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY,
        report_ampc_amr = gene_detection.get_gene_detection_report('ampc', config),
        # Virulence gene detection
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Genomic context investigation
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY,
        # Sequence typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_mlst_mcnally = sequence_typing.get_sequence_typing_report('mlst_mcnally', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_cgmlst_ye = sequence_typing.get_sequence_typing_report('cgmlst_ye', config),
        report_cgmlst_yp = sequence_typing.get_sequence_typing_report('cgmlst_yp', config),
        report_cgmlst_enterobase = sequence_typing.get_sequence_typing_report('cgmlst_enterobase', config),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # Species determination
        report_species = species_determination.OUTPUT_REPORT if 'species' in config['analyses'] else species_determination.OUTPUT_REPORT_EMPTY,
        # Report
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_create_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        input_dict = config['input'],
        input_type = config['input_type'],
        citation_keys = config['citations'],
        detection_method = config['gene_detection']['options']['method']
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
            key_citation=params.citation_keys['main'],
            detection_method=params.detection_method
        ))

        # Set up the report content structure
        report_structure = []
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, input.report_confindr)
        report_structure.extend([
            ('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]),
            ('Species identification', 'species_id', [Path(input.report_rmlst)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_amrfinder, input.report_resfinder4, input.report_ampc_amr)]),
            ('Virulence detection', 'virulence', [Path(input.report_vfdb_core)]),
            ('Genomic context', 'mob_suite', [Path(x) for x in(
                input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_mlst_mcnally, input.report_cgmlst, input.report_cgmlst_ye,
                input.report_cgmlst_yp, input.report_cgmlst_enterobase)]),
            ('Lineage and biotype determination', 'lineage', [Path(input.report_species)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)


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
        quality_checks.OUTPUT_SUMMARY,
        lambda wildcards: contamination_check_kraken.get_summaries(config, wildcards.ext),
        confindr.get_summary(config),
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses'] else [],
        resfinder4.OUTPUT_SUMMARY if 'resfinder4' in config['analyses'] else [],
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='ampc', ext=wildcards.ext) if 'ampc' in config['analyses'] else [],
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # Virulence detection
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='vfdb_core', ext=wildcards.ext) if 'vfdb_core' in config['analyses'] else [],
        # Sequence typing
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst_mcnally', ext=wildcards.ext) if 'mlst_mcnally' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst_ye', ext=wildcards.ext) if 'cgmlst_ye' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst_yp', ext=wildcards.ext) if 'cgmlst_yp' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst_enterobase', ext=wildcards.ext) if 'cgmlst_enterobase' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        species_determination.OUTPUT_SUMMARY if 'species' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.combine_summary_data(input, Path(output.FILE), str(params.ext))

rule link_genomic_context:
    """
    Links the input databases to the genomic context assay.
    """
    input:
        # AMR
        TSV_amrfinder = amrfinder.OUTPUT_TSV if 'amrfinder' in config['analyses'] else [],
        # Virulence
        TSV_gd_vfdb = 'gene_detection/vfdb_core/metadata/tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = gene_detection.OUTPUT_DB_INFORMS.format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))

rule link_species:
    """
    Links the cgST profile matches to the species determination tool.
    """
    input:
        TSV = 'typing/cgmlst/detect_st/tsv.io'
    output:
        TSV_profile_matches = 'species_determination/input/tsv_profile_matches.io',
        TSV_taxonomic = 'species_determination/input/tsv_taxonomic.io'
    params:
        TSV = Path(config['species_determination']['tsv'])
    run:
        import shutil
        from camel.app.snakemake.snakemakeutils import SnakemakeUtils
        from camel.app.io.tooliofile import ToolIOFile
        shutil.copyfile(input.TSV, output.TSV_profile_matches)
        SnakemakeUtils.dump_object([ToolIOFile(params.TSV)], Path(output.TSV_taxonomic))
