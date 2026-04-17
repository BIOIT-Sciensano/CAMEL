from pathlib import Path

from camel.snakefiles import trimming_illumina, trimming_ont, contamination_check_kraken, quality_checks, variant_calling, \
    variant_filtering, gene_detection, sequence_typing, trimming, downsampling, confindr, quast, core, assembly, \
    human_read_scrubbing, read_simulation
from camel.scripts.mycobacteriumpipeline.snakefile import csb_rd, snpit, hsp65, spoligotyping, snplineage, assay51snp, \
    amrdetection
from camel.app.core.snakemake import snakepipelineutils

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
include: confindr.SNAKEFILE
include: contamination_check_kraken.SNAKEFILE
include: quality_checks.SNAKEFILE
include: quast.SNAKEFILE
include: gene_detection.SNAKEFILE
include: sequence_typing.SNAKEFILE
include: variant_calling.SNAKEFILE
include: variant_filtering.SNAKEFILE
include: csb_rd.SNAKEFILE_CSB_RD
include: snpit.SNAKEFILE
include: hsp65.SNAKEFILE
include: assay51snp.SNAKEFILE
include: spoligotyping.SNAKEFILE
include: snplineage.SNAKEFILE
include: amrdetection.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output']['html'],
        config['output']['tsv']

rule report_command_section:
    """
    Creates a HTML report section with the main commands.
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
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config),
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL,
        INFORMS_snpit = snpit.OUTPUT_INFORMS if 'snpit' in config['analyses_selected'] else [],
        INFORMS_16s = str(gene_detection.OUTPUT_INFORMS).format(db='ncbi_16s') if 'ncbi_16s' in config['analyses_selected'] else [],
        INFORMS_csb_rd = str(gene_detection.OUTPUT_INFORMS).format(db='csb_rd') if 'csb_rd' in config['analyses_selected'] else [],
        INFORMS_hsp65 = str(gene_detection.OUTPUT_INFORMS).format(db='hsp65') if 'hsp65' in config['analyses_selected'] else [],
        INFORMS_spoligo = spoligotyping.OUTPUT_INFORMS if 'spoligotyping' in config['analyses_selected'] else [],
        INFORMS_amr = amrdetection.OUTPUT_INFORMS if 'amr' in config['analyses_selected'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

rule report_combine_all:
    """
    Combines the HTML report sections into a single report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = quast.OUTPUT_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = str(quality_checks.OUTPUT_REPORT).format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config),
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_ncbi_16s = gene_detection.get_gene_detection_report('ncbi_16s', config),
        report_51snp = (assay51snp.OUTPUT_REPORT if '51snp' in config['analyses_selected'] else assay51snp.OUTPUT_REPORT_EMPTY),
        report_csb_rd = (csb_rd.OUTPUT_CSB_RD_REPORT if 'csb_rd' in config['analyses_selected'] else csb_rd.OUTPUT_CSB_RD_REPORT_EMPTY),
        report_hsp65 = (hsp65.OUTPUT_REPORT if 'hsp65' in config['analyses_selected'] else hsp65.OUTPUT_REPORT_EMPTY),
        report_snpit = (snpit.OUTPUT_REPORT if 'snpit' in config['analyses_selected'] else snpit.OUTPUT_REPORT_EMPTY),
        # Spoligotyping & lineage determination
        report_spoligo = (spoligotyping.OUTPUT_REPORT if 'spoligotyping' in config['analyses_selected'] else spoligotyping.OUTPUT_REPORT_EMPTY),
        report_snp_lineage = (snplineage.OUTPUT_REPORT if 'snp_lineage' in config['analyses_selected'] else snplineage.OUTPUT_REPORT_EMPTY),
        # AMR
        report_amr = (amrdetection.OUTPUT_REPORT if 'amr' in config['analyses_selected'] else amrdetection.OUTPUT_REPORT_EMPTY),
        report_amr_genes = amrdetection.OUTPUT_REPORT_CDS,
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        # Report
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        input_dict = config['input'],
        citation_keys = config['citations'],
        gene_detection_method = config['gene_detection']['options']['method'],
        typing_method = config['sequence_typing']['options']['method']
    run:
        import datetime
        from camel.app.scriptutils.basepipe import basepipeutils
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput
        from camel.app.scriptutils.model import InputType

        # Add the header section
        script_input: ScriptInput = ScriptInput.from_dict(params.input_dict)
        report = snakepipelineutils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        section = snakepipelineutils.create_input_section(
            sample_name=script_input.name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=script_input.input_str,
            input_type=script_input.type_.value,
            key_citation=params.citation_keys['main'],
            extra_data=[
                ('Gene detection method', params.gene_detection_method),
                ('Typing method', params.typing_method),
            ]
        )
        if script_input.type_ == InputType.FASTA:
            section.add_warning_message(
                'SNP-based assays are run on simulated reads from the assembled contigs, which may differ from the '
                'original reads.')
        report.add_html_object(section)

        # Add report content
        report_structure = []
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))

        # Add output sections
        report_structure.extend([
            ('Variant calling', 'variant', [Path(input.report_variant)]),
            ('Species identification', 'identification', [
                Path(input.report_rmlst), Path(input.report_ncbi_16s), Path(input.report_snpit),
                Path(input.report_csb_rd), Path(input.report_hsp65), Path(input.report_51snp)]),
            ('Spoligotyping and lineage', 'spoligotyping', [
                Path(input.report_spoligo), Path(input.report_snp_lineage)]),
            ('AMR detection', 'amr', [Path(input.report_amr), Path(input.report_amr_genes)]),
            ('Sequence typing', 'typing', [Path(input.report_mlst), Path(input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        snakepipelineutils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    Combines the summary output files of the different assays into a single summary file.
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
        variant_calling.get_summaries(config),
        variant_filtering.OUTPUT_SUMMARY,
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='ncbi_16s', ext=wildcards.ext) if 'ncbi_16s' in config['analyses_selected'] else [],
        csb_rd.OUTPUT_CSB_RD_SUMMARY if 'csb_rd' in config['analyses_selected'] else [],
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='hsp65', ext=wildcards.ext) if 'hsp65' in config['analyses_selected'] else [],
        assay51snp.OUTPUT_SUMMARY if '51snp' in config['analyses_selected'] else [],
        snpit.OUTPUT_SUMMARY if 'snpit' in config['analyses_selected'] else [],
        spoligotyping.OUTPUT_SUMMARY if 'spoligotyping' in config['analyses_selected'] else [],
        snplineage.OUTPUT_SUMMARY if 'snp_lineage' in config['analyses_selected'] else [],
        amrdetection.OUTPUT_SUMMARY if 'amr' in config['analyses_selected'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses_selected'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses_selected'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses_selected'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.combine_summary_data(input, Path(output.FILE), str(params.ext))
