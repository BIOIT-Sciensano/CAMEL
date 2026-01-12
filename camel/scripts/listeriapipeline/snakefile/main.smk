from pathlib import Path

from camel.snakefiles import trimming_illumina, gene_detection, trimming, contamination_check_kraken, \
    quality_checks, variant_calling, variant_filtering, sequence_typing, downsampling, confindr, quast, amrfinder, \
    core, trimming_ont, assembly, human_read_scrubbing, resfinder4, mobsuite, read_simulation

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
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: mobsuite.SNAKEFILE
include: gene_detection.SNAKEFILE
include: sequence_typing.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    Rule to generate the required output files.
    """
    input:
        HTML = config['output']['html'],
        TSV = config['output']['tsv'],
        JSON = config['output']['json'] if config['output'].get('json') is not None else []

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
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input']['type']),
        # Variant calling
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        # AMRFinder
        INFORMS_amrfnder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        # Gene detection
        INFORMS_resfinder = gene_detection.OUTPUT_INFORMS.format(db='resfinder') if 'resfinder' in config['analyses'] else [],
        INFORMS_virulencefinder = gene_detection.OUTPUT_INFORMS.format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_vfdb_core = gene_detection.OUTPUT_INFORMS.format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        INFORMS_plasmidfinder = gene_detection.OUTPUT_INFORMS.format(db='plasmidfinder') if 'plasmidfinder' in config['analyses'] else [],
        # Sequence typing
        INFORMS_rmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        INFORMS_mlst = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_species = sequence_typing.OUTPUT_INFORMS.format(scheme='species_confirmation') if 'species_confirmation' in config['analyses'] else [],
        INFORMS_cgmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_typing_amr = sequence_typing.OUTPUT_INFORMS.format(scheme='typing_amr') if 'typing_amr' in config['analyses'] else [],
        INFORMS_virulence = sequence_typing.OUTPUT_INFORMS.format(scheme='typing_virulence') if 'typing_virulence' in config['analyses'] else [],
        INFORMS_pcr_sero = sequence_typing.OUTPUT_INFORMS.format(scheme='pcr_serogroup') if 'pcr_serogroup' in config['analyses'] else [],
        # Plasmid characterization
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS if 'mob_suite' in config['analyses'] else[]
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

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
        report_adv_qc = str(quality_checks.OUTPUT_REPORT).format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        # AMR detection
        report_amrfinder = (amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY),
        report_resfinder4 = (resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY),
        # Virulence detection
        report_virulence = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY,
        # Typing
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_species = sequence_typing.get_sequence_typing_report('species_confirmation', config),
        report_amr_typing = sequence_typing.get_sequence_typing_report('typing_amr', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_pcr_serogroup = sequence_typing.get_sequence_typing_report('pcr_serogroup', config),
        report_viru_typing = sequence_typing.get_sequence_typing_report('typing_virulence', config),
        # Report
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        input_dict = config['input'],
        gene_detection_method = config['gene_detection']['options']['method'],
        typing_method = config['sequence_typing']['options']['method']
    run:
        import datetime
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput
        from camel.app.scriptutils.basepipe import basepipeutils
        from camel.app.core.snakemake import snakepipelineutils

        # Add the header section
        script_input = ScriptInput.from_dict(params.input_dict)
        report = snakepipelineutils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(snakepipelineutils.create_input_section(
            sample_name=script_input.name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=script_input.input_str,
            input_type=script_input.type_.value,
            extra_data=[
                ('Gene detection method', params.gene_detection_method),
                ('Typing method', params.typing_method),
            ]
        ))

        # Set up the report content structure
        report_structure = []
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value,input.reports_downsampling,input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination,input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        report_structure.extend([
            ('Species identification', 'species', [Path(x) for x in (
                input.report_rmlst, input.report_species, input.report_pcr_serogroup)]),
            ('AMR detection', 'amr', [Path(x) for x in (input.report_amrfinder, input.report_resfinder4, input.report_amr_typing)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_virulence, input.report_vfdb_core, input.report_viru_typing)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'typing', [Path(x) for x in (input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        snakepipelineutils.add_report_content(report, report_structure)

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
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='virulencefinder', ext=wildcards.ext) if 'virulencefinder' in config['analyses'] else [],
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='vfdb_core', ext=wildcards.ext) if 'vfdb_core' in config['analyses'] else [],
        # Plasmid characterization
        lambda wildcards: str(gene_detection.OUTPUT_SUMMARY).format(db='plasmidfinder', ext=wildcards.ext) if 'plasmidfinder' in config['analyses'] else [],
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # Sequence typing
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='species_confirmation', ext=wildcards.ext) if 'species_confirmation' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='typing_amr', ext=wildcards.ext) if 'typing_amr' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='typing_virulence', ext=wildcards.ext) if 'typing_amr' in config['analyses'] else [],
        lambda wildcards: str(sequence_typing.OUTPUT_SUMMARY).format(scheme='pcr_serogroup', ext=wildcards.ext) if 'pcr_serogroup' in config['analyses'] else [],
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
        TSV_gd_virulencefinder = 'gene_detection/virulencefinder/metadata/tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_virulencefinder = str(gene_detection.OUTPUT_DB_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
