from pathlib import Path

from camel.snakefiles import trimming_illumina, trimming_ont, gene_detection, trimming, contamination_check_kraken, \
    quality_checks, variant_calling, variant_filtering, sequence_typing, lrefinder, downsampling, confindr, quast, \
    core, assembly, amrfinder, resfinder4, bacmet, mobsuite, human_read_scrubbing, read_simulation
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
include: quast.SNAKEFILE
include: contamination_check_kraken.SNAKEFILE
include: confindr.SNAKEFILE
include: quality_checks.SNAKEFILE
include: variant_calling.SNAKEFILE
include: variant_filtering.SNAKEFILE
include: lrefinder.SNAKEFILE
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: gene_detection.SNAKEFILE
include: mobsuite.SNAKEFILE
include: bacmet.SNAKEFILE
include: sequence_typing.SNAKEFILE


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
        INFORMS_simulation = read_simulation.OUTPUT_INFORMS if config['input']['type'] == 'fasta' else [],
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = quast.OUTPUT_INFORMS,
        INFORMS_busco = quast.OUTPUT_INFORMS_BUSCO,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input']['type']),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        INFORMS_lrefinder = lrefinder.OUTPUT_INFORMS if 'lrefinder' in config['analyses'] else [],
        INFORMS_amrfinder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS if 'mob_suite' in config['analyses'] else [],
        INFORMS_virulencefinder = gene_detection.OUTPUT_INFORMS.format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        INFORMS_vfdb_core = gene_detection.OUTPUT_INFORMS.format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        IFNORMS_bacmet = bacmet.OUTPUT_INFORMS if 'bacmet' in config['analyses'] else [],
        INFORMS_prodigal = bacmet.OUTPUT_PRODIGAL_INFORMS if 'bacmet' in config['analyses'] else [],
        INFORMS_mlst = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst') if 'mlst' in config['analyses'] else []
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
        report_adv_qc = quality_checks.OUTPUT_REPORT.format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        # Species identification
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # AMR detection
        report_lrefinder = lrefinder.OUTPUT_REPORT if 'lrefinder' in config['analyses'] else lrefinder.OUTPUT_REPORT_EMPTY,
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_resfinder4 = resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY,
        # Virulence gene detection
        report_virulencefinder = gene_detection.get_gene_detection_report('virulencefinder', config),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        # Plasmid characterization
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY,
        # BacMet
        report_prodigal = bacmet.OUTPUT_PRODIGAL_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_PRODIGAL_REPORT_EMPTY,
        report_bacmet = bacmet.OUTPUT_REPORT if 'bacmet' in config['analyses'] else bacmet.OUTPUT_REPORT_EMPTY,
        # Typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config) if not config['is_generic'] else [],
        report_mlst_bezdicek = sequence_typing.get_sequence_typing_report('mlst_bezdicek', config) if not config['is_generic'] else [],
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config) if not config['is_generic'] else [],
        # Report
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        species = config['selected_species'],
        input_dict = config['input'],
        citation_keys = config['citations'],
        detection_method = config['gene_detection']['options']['method']
    run:
        import datetime
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput
        from camel.app.scriptutils.basepipe import basepipeutils

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
            detection_method=params.detection_method,
            extra_data=[('Selected species', f'<i>{params.species}</i>')],
        ))

        # Add report content
        report_structure = []
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value,input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))

        # Typing (additional MLST scheme for E. faecium)
        if params.species == 'Enterococcus faecalis':
            reports_typing = (Path(input.report_mlst), Path(input.report_cgmlst))
        elif params.species == 'Enterococcus faecium':
            reports_typing = (Path(input.report_mlst), Path(input.report_mlst_bezdicek), Path(input.report_cgmlst))
        else:
            reports_typing = ()

        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('AMR detection', 'amr', [Path(x) for x in (
                input.report_lrefinder, input.report_amrfinder, input.report_resfinder4)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (
                input.report_virulencefinder, input.report_vfdb_core)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Biocide and metal resistance', 'bacmet', [Path(input.report_prodigal), Path(input.report_bacmet)])
        ])
        if len(reports_typing) > 0:
            report_structure.append(('Sequence typing', 'typing', reports_typing))
        report_structure.extend([
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
        lambda wildcards: downsampling.get_summaries(config, ext=wildcards.ext),
        trimming.get_summaries(config),
        quast.OUTPUT_SUMMARY,
        lambda wildcards: contamination_check_kraken.get_summaries(config, ext=wildcards.ext),
        confindr.get_summary(config),
        quality_checks.OUTPUT_SUMMARY,
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses'] else [],
        # AMR detection
        lrefinder.OUTPUT_SUMMARY if 'lrefinder' in config['analyses'] else [],
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses'] else [],
        resfinder4.OUTPUT_SUMMARY if 'resfinder4' in config['analyses'] else [],
        # Virulence detection
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='virulencefinder', ext=wildcards.ext) if 'virulencefinder' in config['analyses'] else [],
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='vfdb_core', ext=wildcards.ext) if 'vfdb_core' in config['analyses'] else [],
        # Plasmid characterization
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='plasmidfinder', ext=wildcards.ext) if 'plasmidfinder' in config['analyses'] else [],
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses'] else [],
        # BacMet
        bacmet.OUTPUT_SUMMARY if 'bacmet' in config['analyses'] else [],
        # Sequence typing
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst_bezdicek', ext=wildcards.ext) if 'mlst_bezdicek' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses'] else []
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
        TSV_gd_vfdb = 'gene_detection/vfdb_core/metadata/tsv.io' if 'vfdb_core' in config['analyses'] else [],
        INFORMS_gd_vfdb = str(gene_detection.OUTPUT_DB_INFORMS).format(db='vfdb_core') if 'vfdb_core' in config['analyses'] else [],
        TSV_gd_virulencefinder = 'gene_detection/virulencefinder/metadata/tsv.io' if 'virulencefinder' in config['analyses'] else [],
        INFORMS_gd_virulencefinder = str(gene_detection.OUTPUT_DB_INFORMS).format(db='virulencefinder') if 'virulencefinder' in config['analyses'] else [],
        # BacMet
        TSV_bacmet = 'bacmet/hit_filtering/tsv.io' if 'bacmet' in config['analyses'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
