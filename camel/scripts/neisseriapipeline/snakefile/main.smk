from pathlib import Path

from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.snakefiles import trimming, trimming_illumina, quality_checks, variant_calling, variant_filtering, \
    contamination_check_kraken, sequence_typing, downsampling, confindr, quast, core, trimming_ont, \
    assembly, human_read_scrubbing, amrfinder, resfinder4, read_simulation
from camel.scripts.neisseriapipeline.snakefile import serogroup_determination, gmats, mendevar


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
include: amrfinder.SNAKEFILE
include: resfinder4.SNAKEFILE
include: gmats.SNAKEFILE
include: mendevar.SNAKEFILE
include: serogroup_determination.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        config['output']['html'],
        config['output']['tsv'],
        config['output']['json'] if config['output'].get('json') is not None else []

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
        INFORMS_assembly_map = assembly.get_qc_informs(config, config['input']['type']),
        INFORMS_variant_calling_all = variant_calling.get_command_informs(config) if 'variant_calling' in config['analyses'] else [],
        INFORMS_variant_filtering_all = variant_filtering.OUTPUT_INFORMS_ALL if 'variant_calling' in config['analyses'] else [],
        INFORMS_amrfnder = amrfinder.OUTPUT_INFORMS if 'amrfinder' in config['analyses'] else [],
        INFORMS_resfinder4 = resfinder4.OUTPUT_INFORMS if 'resfinder4' in config['analyses'] else [],
        INFORMS_rmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='rmlst') if 'rmlst' in config['analyses'] else [],
        INFORMS_mlst = sequence_typing.OUTPUT_INFORMS.format(scheme='mlst') if 'mlst' in config['analyses'] else [],
        INFORMS_rplf = sequence_typing.OUTPUT_INFORMS.format(scheme='rplf') if 'rplf' in config['analyses'] else [],
        INFORMS_bast = sequence_typing.OUTPUT_INFORMS.format(scheme='bast') if 'bast' in config['analyses'] else [],
        INFORMS_pora = sequence_typing.OUTPUT_INFORMS.format(scheme='pora') if 'pora' in config['analyses'] else [],
        INFORMS_porb = sequence_typing.OUTPUT_INFORMS.format(scheme='porb') if 'porb' in config['analyses'] else [],
        INFORMS_feta = sequence_typing.OUTPUT_INFORMS.format(scheme='feta') if 'feta' in config['analyses'] else [],
        INFORMS_amr = sequence_typing.OUTPUT_INFORMS.format(scheme='resistance_genes') if 'resistance_genes' in config['analyses'] else [],
        INFORMS_vaccine = sequence_typing.OUTPUT_INFORMS.format(scheme='vaccine_targets') if 'vaccine_targets' in config['analyses'] else [],
        INFORMS_fhbp = sequence_typing.OUTPUT_INFORMS.format(scheme='fhbp') if 'fhbp' in config['analyses'] else [],
        INFORMS_cgmlst = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst') if 'cgmlst' in config['analyses'] else [],
        INFORMS_cgmlst_v3 = sequence_typing.OUTPUT_INFORMS.format(scheme='cgmlst_v3') if 'cgmlst_v3' in config['analyses'] else [],
        INFORMS_serogroup = serogroup_determination.OUTPUT_INFORMS if 'serogroup' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

rule neisseria_additional_resistance_gene_metadata:
    """
    This rule is used to add resistance gene metadata for penA and rpoB genes.
    The data is parsed from the PubMLST webpage.
    """
    input:
        hits = sequence_typing.OUTPUT_HITS.format(
            scheme='resistance_genes',
            locus_type='DNA',
            detection_method=config['sequence_typing']['options']['method']),
        VAL_HTML = sequence_typing.get_sequence_typing_report('resistance_genes', config),
        INFORMS_scheme = sequence_typing.OUTPUT_DB_INFORMS.format(scheme='resistance_genes')
    output:
        VAL_HTML = 'typing/resistance_genes/metadata/html.iob'
    params:
        dir_ = 'typing/resistance_genes/metadata',
        loci='penA, rpoB'
    run:
        from camel.app.core.snakemake.step import Step
        from camel.app.tools.pipelines.neisseria.resistancemetadataextractor import ResistanceMetadataExtractor
        extractor = ResistanceMetadataExtractor()
        snakemakeutils.add_pickle_inputs(extractor, input)
        step = Step(rule_name=str(rule), tool=extractor, dir_=Path(str(params.dir_)))
        extractor.update_parameters(loci=params.loci)
        step.run()
        snakemakeutils.dump_tool_outputs(extractor, output)

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
        report_adv_qc= quality_checks.OUTPUT_REPORT.format(input_type=config['input']['type']),
        report_variant = variant_calling.get_reports(config) if 'variant_calling' in config['analyses'] else [],
        report_rmlst = sequence_typing.get_sequence_typing_report('rmlst', config),
        # AMR detection
        report_amrfinder = amrfinder.OUTPUT_REPORT if 'amrfinder' in config['analyses'] else amrfinder.OUTPUT_REPORT_EMPTY,
        report_resfinder4 = resfinder4.OUTPUT_REPORT if 'resfinder4' in config['analyses'] else resfinder4.OUTPUT_REPORT_EMPTY,
        # Sequence typing
        report_mlst = sequence_typing.get_sequence_typing_report('mlst', config),
        report_cgmlst = sequence_typing.get_sequence_typing_report('cgmlst', config),
        report_cgmlst_v3 = sequence_typing.get_sequence_typing_report('cgmlst_v3', config),
        report_pora = sequence_typing.get_sequence_typing_report('pora', config),
        report_porb = sequence_typing.get_sequence_typing_report('porb', config),
        report_feta = sequence_typing.get_sequence_typing_report('feta', config),
        report_rplf = sequence_typing.get_sequence_typing_report('rplf', config),
        report_vaccine_targets = sequence_typing.get_sequence_typing_report('vaccine_targets', config),
        report_resistance_genes = rules.neisseria_additional_resistance_gene_metadata.output.VAL_HTML if 'resistance_genes' in config['analyses'] else sequence_typing.OUTPUT_REPORT_EMPTY.format(scheme='resistance_genes'),
        report_fhbp = sequence_typing.get_sequence_typing_report('fhbp', config),
        report_bast = sequence_typing.get_sequence_typing_report('bast', config),
        report_gmats = gmats.OUTPUT_REPORT if 'gmats' in config['analyses'] else gmats.OUTPUT_REPORT_EMPTY,
        report_mendevar = mendevar.OUTPUT_REPORT if 'mendevar' in config['analyses'] else mendevar.OUTPUT_REPORT_EMPTY,
        report_serogroup = serogroup_determination.OUTPUT_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_REPORT_EMPTY,
        report_serogroup_legacy = serogroup_determination.OUTPUT_LEGACY_REPORT if 'serogroup' in config['analyses'] else serogroup_determination.OUTPUT_LEGACY_REPORT_EMPTY,
        report_citations = core.OUTPUT_HTML_CITATIONS,
        report_commands = rules.report_create_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        input_dict = config['input'],
        typing_method = config['sequence_typing']['options']['method'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.scriptutils.basepipe import basepipeutils
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
            extra_data=[('Typing method', params.typing_method)],
            key_citation=params.citation_keys['main'],
        ))

        # Set up the report content structure
        report_structure = []
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        basepipeutils.add_content_trim_basic_qc(
            report_structure, script_input.type_.value, input.reports_downsampling, input.reports_trimming)
        report_structure.append(('Assembly', 'assembly', [Path(input.report_quast)]))
        basepipeutils.add_content_contamination_check(
            report_structure, script_input.type_.value, input.reports_contamination, input.report_confindr)
        report_structure.append(('Advanced QC', 'adv_qc', [Path(input.report_adv_qc)]))
        if 'variant_calling' in config['analyses']:
            report_structure.append(('Variant calling', 'variant', [Path(input.report_variant)]))
        report_structure.extend([
            ('Species identification', 'species', [Path(input.report_rmlst)]),
            ('Antimicrobial resistance characterization', 'amr', [Path(x) for x in (
                input.report_amrfinder, input.report_resfinder4)]),
            ('Sequence typing', 'st', [Path(x) for x in (
                input.report_mlst, input.report_rplf, input.report_pora, input.report_porb, input.report_feta,
                input.report_resistance_genes, input.report_vaccine_targets, input.report_fhbp, input.report_cgmlst,
                input.report_cgmlst_v3)]),
            ('Antigen typing', 'at', [Path(x) for x in (input.report_bast, input.report_gmats, input.report_mendevar)]),
            ('Serogroup determination', 'serogroup', [Path(
                input.report_serogroup), Path(input.report_serogroup_legacy)]),
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
        quality_checks.OUTPUT_SUMMARY,
        lambda wildcards: contamination_check_kraken.get_summaries(config, wildcards.ext),
        confindr.get_summary(config),
        variant_calling.get_summaries(config) if 'variant_calling' in config['analyses'] else [],
        amrfinder.OUTPUT_SUMMARY if 'amrfinder' in config['analyses'] else [],
        resfinder4.OUTPUT_SUMMARY if 'resfinder4' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='rmlst', ext=wildcards.ext) if 'rmlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='mlst', ext=wildcards.ext) if 'mlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='rplf', ext=wildcards.ext) if 'rplf' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='bast', ext=wildcards.ext) if 'bast' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='pora', ext=wildcards.ext) if 'pora' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='porb', ext=wildcards.ext) if 'porb' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='feta', ext=wildcards.ext) if 'feta' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='fhbp', ext=wildcards.ext) if 'fhbp' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='resistance_genes', ext=wildcards.ext) if 'resistance_genes' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='vaccine_targets', ext=wildcards.ext) if 'vaccine_targets' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='cgmlst', ext=wildcards.ext) if 'cgmlst' in config['analyses'] else [],
        lambda wildcards: sequence_typing.OUTPUT_SUMMARY.format(scheme='cgmlst_v3', ext=wildcards.ext) if 'cgmlst_v3' in config['analyses'] else [],
        gmats.OUTPUT_SUMMARY if 'gmats' in config['analyses'] else [],
        mendevar.OUTPUT_SUMMARY if 'mendevar' in config['analyses'] else [],
        serogroup_determination.OUTPUT_SUMMARY if 'serogroup' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.combine_summary_data(input, Path(output.FILE), str(params.ext))
