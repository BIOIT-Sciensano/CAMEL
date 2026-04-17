from pathlib import Path

from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.snakefiles import core, assembly, downsampling, quast, confindr, trimming, trimming_illumina, \
    quality_checks, variant_calling, variant_filtering, contamination_check_kraken, sequence_typing, amrfinder, \
    trimming_ont, gene_detection, mobsuite, human_read_scrubbing, read_simulation
from camel.scripts.bacilluspipeline.snakefile import btyper, ani, straingst

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
include: sequence_typing.SNAKEFILE
include: btyper.SNAKEFILE
include: amrfinder.SNAKEFILE
include: gene_detection.SNAKEFILE
include: mobsuite.SNAKEFILE
include: ani.SNAKEFILE
include: straingst.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        HTML = config['output']['html'],
        TSV = config['output']['tsv']

#####################################
# Linking workflow inputs & outputs #
#####################################
rule link_fasta_to_tools_subtilis:
    """
    This rule links the output of the assembly workflow to the fastANI workflow if the species is B. subtilis.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = ani.INPUT_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA}"

rule link_fasta_to_tools_cereus:
    """
    This rule links the output of the assembly workflow to the BTyper workflow if the species is B. cereus.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = btyper.INPUT_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA}"

rule main_update_gmm_report:
    """
    This rule updates the report of the GMM detection assay with an interpretation paragraph.
    """
    input:
        VAL_HTML_VECTORS = gene_detection.get_gene_detection_report('gmm_genes_vectors', config, analysis_name='gmo'),
        VAL_HTML_JUNCTIONS = gene_detection.get_gene_detection_report('gmm_junctions', config, analysis_name='gmo'),
        TSV_STRAINS = straingst.get_summaries(config, ext='tsv'),
        TSV_GMM_VECTORS = gene_detection.OUTPUT_SUMMARY.format(db='gmm_genes_vectors', ext='tsv') if 'gmo' in config['analyses_selected'] else [],
        TSV_GMM_JUNCTIONS = gene_detection.OUTPUT_SUMMARY.format(db='gmm_junctions', ext='tsv') if 'gmo' in config['analyses_selected'] else [],
        TSV_GMM_DB = config['gene_detection']['dbs']['gmm_genes_vectors']['known_gmm_constructs']
    output:
        VAL_HTML = 'gene_detection/gmo/updated_html_report.iob'
    params:
        running_dir = 'gene_detection/gmo'
    run:
        from camel.app.tools.pipelines.bacillus.updategmmreport import UpdateGMMReport
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.core.snakemake.step import Step
        gmmupdater = UpdateGMMReport()
        snakemakeutils.add_io_inputs(gmmupdater, input, excluded_keys=['TSV_STRAINS', 'TSV_GMM_VECTORS', 'TSV_GMM_JUNCTIONS', 'TSV_GMM_DB'])
        gmmupdater.add_input_files({
            'TSV_STRAINS': [ToolIOFile(Path(x)) for x in input.TSV_STRAINS],
            'TSV_GMM_DB': [ToolIOFile(Path(input.TSV_GMM_DB))],
            'TSV_GMM_VECTORS': [ToolIOFile(Path(input.TSV_GMM_VECTORS))],
            'TSV_GMM_JUNCTIONS': [ToolIOFile(Path(input.TSV_GMM_JUNCTIONS))]
             })
        step = Step(rule_name=str(rule), tool=gmmupdater, dir_=Path(str(params.running_dir)))
        step.run()
        snakemakeutils.dump_io_outputs(gmmupdater, output)

##########
# Report #
##########
rule report_create_commands_section:
    """
    Creates the section with the commands.
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
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses_selected'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses_selected'] else [],
        INFORMS_assembly_map = assembly.get_qc_informs(config['input']['type']),
        INFORMS_btyper = btyper.OUTPUT_INFORMS if 'btyper' in config['analyses_selected'] else [],
        INFORMS_fastani = ani.OUTPUT_INFORMS if 'fastani' in config['analyses_selected'] else [],
        INFORMS_straingst = straingst.get_command_informs(config),
        INFORMS_amrfinder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses_selected'] else [],
        INFORMS_vfdb_core = gene_detection.OUTPUT_INFORMS.format(db='vfdb_core') if 'vfdb_core' in config['analyses_selected'] else [],
        INFORMS_plasmidfinder = gene_detection.OUTPUT_INFORMS.format(db='plasmidfinder') if 'plasmidfinder' in config['analyses_selected'] else [],
        INFORMS_mob_suite = mobsuite.OUTPUT_INFORMS if 'mob_suite' in config['analyses_selected'] else [],
        INFORMS_mlst_cereus = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst_cereus') if 'mlst_cereus' in config['analyses_selected'] else [],
        INFORMS_mlst_subtilis = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst_subtilis') if 'mlst_subtilis' in config['analyses_selected'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

rule report_content_cereus:
    """
    Creates the main content of the report when the detected species is Bacillus cereus.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = quast.OUTPUT_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = quality_checks.OUTPUT_REPORT.format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses_selected'] else [],
        report_btyper = btyper.OUTPUT_REPORT if 'btyper' in config['analyses_selected'] else btyper.OUTPUT_REPORT_EMPTY,
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses_selected'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mob_suite' in config['analyses_selected'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses_selected'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY,
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst_cereus', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst_cereus', config),
        report_commands = rules.report_create_commands_section.output.HTML,
        report_citations = core.OUTPUT_HTML_CITATIONS
    output:
        HTML = 'report/report_cereus.html'
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        species = config['species'][config['species_selected']]['full_name'],
        input_dict = config['input'],
        citation_keys = config['citations'],
        gene_detection_method = config['gene_detection']['options']['method'],
        typing_method = config['sequence_typing']['options']['method']
    run:
        import datetime
        from camel.app.scriptutils.basepipe import basepipeutils
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput

        # Add the header section+
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
                ('Selected species', f'<i>{params.species}</i>'),
                ('Gene detection method', params.gene_detection_method),
                ('Typing method', params.typing_method),
            ],
            key_citation=params.citation_keys['main'],
        ))
        report_structure = []

        # Core sections (shared)
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses_selected']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        # Custom assays (B. cereus)
        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('BTyper3', 'btyper3', [Path(input.report_btyper)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_mlst, input.report_cgmlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        snakepipelineutils.add_report_content(report, report_structure)
        report.save()

rule report_content_subtilis:
    """
    Creates the main content of the report when the detected species is Bacillus subtilis.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = quast.OUTPUT_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = quality_checks.OUTPUT_REPORT.format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses_selected'] else [],
        report_fastani = ani.OUTPUT_REPORT if 'fastani' in config['analyses_selected'] else ani.OUTPUT_REPORT_EMPTY,
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses_selected'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_gmo = rules.main_update_gmm_report.output.VAL_HTML if 'gmo' in config['analyses_selected'] else gene_detection.get_gene_detection_report('gmm_genes_vectors', config, analysis_name='gmo'),
        report_junctions = gene_detection.get_gene_detection_report('gmm_junctions', config, analysis_name='gmo'),
        report_vfdb_core = gene_detection.get_gene_detection_report('vfdb_core', config),
        report_plasmidfinder = gene_detection.get_gene_detection_report('plasmidfinder', config),
        report_mob_suite = mobsuite.OUTPUT_REPORT if 'mobsuite' in config['analyses_selected'] else mobsuite.OUTPUT_REPORT_EMPTY,
        report_genomic_context = mobsuite.OUTPUT_CONTEXT_REPORT if 'mob_suite' in config['analyses_selected'] else mobsuite.OUTPUT_CONTEXT_REPORT_EMPTY,
        reports_straingst = straingst.get_reports(config),
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        report_mlst = sequence_typing.get_sequence_typing_report('mlst_subtilis', config),
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_create_commands_section.output.HTML
    output:
        HTML = 'report/report_subtilis.html'
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        species = config['species'][config['species_selected']]['full_name'],
        input_dict = config['input'],
        citation_keys = config['citations'],
        gene_detection_method = config['gene_detection']['options']['method'],
        typing_method = config['sequence_typing']['options']['method']
    run:
        import datetime
        from camel.app.scriptutils.basepipe import basepipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput

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
                ('Selected species', f'<i>{params.species}</i>'),
                ('Gene detection method', params.gene_detection_method),
                ('Typing method', params.typing_method),
            ],
            key_citation=params.citation_keys['main'],
        ))

        # Create the report
        report_structure = []

        # Core sections
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses_selected']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        # B. subtilis assays
        report_structure.extend([
            ('FastANI', 'fastani', [Path(input.report_fastani)]),
            ('StrainGST', 'straingst', [Path(x) for x in input.reports_straingst]),
            ('GMO detection', 'gmo', [Path(input.report_junctions), Path(input.report_gmo)]),
            # ('GMO detection', 'gmo', [Path(x) for x in (input.report_gmo,)]),
            ('Virulence detection', 'virulence', [Path(x) for x in (input.report_vfdb_core,)]),
            ('AMRFinder results', 'amrfinder', [Path(input.report_amrfinder)]),
            ('Plasmid characterization', 'plasmid', [Path(x) for x in (
                input.report_plasmidfinder, input.report_mob_suite, input.report_genomic_context)]),
            ('Sequence typing', 'st', [Path(x) for x in (input.report_rmlst, input.report_mlst)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ])
        snakepipelineutils.add_report_content(report, report_structure)
        report.save()

rule report_select_by_species:
    """
    Selects the report content based on the detected species.
    """
    input:
        HTML = f'report/report_{config["species_selected"]}.html'
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir']
    shell:
        """
        cp {input.HTML} {output.HTML}
        """

rule summary_combine_all:
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
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses_selected'] else [],
        lambda wildcards: straingst.get_summaries(config, ext=wildcards.ext),
        btyper.OUTPUT_SUMMARY if 'btyper' in config['analyses_selected'] else [],
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses_selected'] else [],
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='vfdb_core', ext=wildcards.ext) if 'vfdb_core' in config['analyses_selected'] else [],
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='gmm_genes_vectors', ext=wildcards.ext) if 'gmo' in config['analyses_selected'] else [],
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='gmm_junctions', ext=wildcards.ext) if 'gmo' in config['analyses_selected'] else [],
        lambda wildcards: gene_detection.OUTPUT_SUMMARY.format(db='plasmidfinder', ext=wildcards.ext) if 'plasmidfinder' in config['analyses_selected'] else [],
        ani.OUTPUT_SUMMARY if 'fastani' in config['analyses_selected'] else [],
        mobsuite.OUTPUT_SUMMARY if 'mob_suite' in config['analyses_selected'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses_selected'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst_cereus', ext=wildcards.ext) if 'mlst_cereus' in config['analyses_selected'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst_subtilis', ext=wildcards.ext) if 'mlst_subtilis' in config['analyses_selected'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='cgmlst_cereus', ext=wildcards.ext) if 'cgmlst_cereus' in config['analyses_selected'] else []
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
        TSV_amrfinder = amrfinder.OUTPUT_TSV if 'amrfinder' in config['analyses_selected'] else [],
        # Virulence
        TSV_gd_vfdb = 'gene_detection/vfdb_core/metadata/tsv.io' if 'vfdb_core' in config['analyses_selected'] else [],
        INFORMS_gd_vfdb = gene_detection.OUTPUT_DB_INFORMS.format(db='vfdb_core') if 'vfdb_core' in config['analyses_selected'] else []
    output:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io'
    run:
        mobsuite.collect_genomic_context_input(input, Path(output.TSV), Path(output.INFORMS))
