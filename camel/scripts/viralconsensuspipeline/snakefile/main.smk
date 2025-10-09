import shutil
from pathlib import Path

from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_illumina, trimming_ont, trimming, downsampling, \
    contamination_check_kraken, core, human_read_scrubbing, assembly
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping, refselection, preprocess, \
    multiallelicsites, nextclade3, antivirals

#######################
# Included snakefiles #
#######################
include: core.SNAKEFILE
include: assembly.SNAKEFILE
include: human_read_scrubbing.SNAKEFILE
include: downsampling.SNAKEFILE
include: trimming_ont.SNAKEFILE
include: trimming_illumina.SNAKEFILE
include: contamination_check_kraken.SNAKEFILE
include: refselection.SNAKEFILE
include: preprocess.SNAKEFILE
include: iterativemapping.SNAKEFILE
include: nextclade3.SNAKEFILE_NEXTCLADE
include: multiallelicsites.SNAKEFILE
include: antivirals.SNAKEFILE

#########
# Rules #
#########
rule all:
    """
    This rules ensures that the required output files are generated.
    """
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular'],
        JSON = config['output_json'] if config['output_json'] is not None else []

rule link_fasta_to_iterative_mapping:
    """
    Selects the FASTA file used as a reference for read mapping.
    """
    input:
        FASTA = refselection.OUTPUT_FASTA if config['fasta_ref'] is None else []
    output:
        FASTA = iterativemapping.INPUT_FASTA_REF
    params:
        input_type = config['input_type'],
        fasta_ref = config['fasta_ref']
    run:
        from camel.app.io.tooliofile import ToolIOFile
        if params.fasta_ref is not None:
            snakemakeutils.dump_object([ToolIOFile(Path(params.fasta_ref))], Path(output.FASTA))
        else:
            shutil.copyfile(input.FASTA, output.FASTA)

rule select_fasta_file:
    """
    This rule selects the fasta file to send to other workflows.
    """
    input:
        FASTA = iterativemapping.get_fasta(config)
    output:
        FASTA = 'fasta.io'
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule report_create_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config),
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_reference_selection = Path('ref_selection') / 'mash_screen' / refselection.get_segments(
            Path(config['ref_selection']['db']))[0] / 'informs.io' if config['fasta_ref'] is None and config['input_type'] != 'fasta' and config['ref_selection'].get('db') is not None else [],
        INFORMS_preprocess = preprocess.OUTPUT_INFORMS if config['input_type'] != 'fasta' else [],
        INFORMS_iterative_mapping = iterativemapping.OUTPUT_INFORMS if config['input_type'] != 'fasta' else [],
        INFORMS_mash = nextclade3.OUTPUT_INFORMS_MASH if config['nextclade'].get('db_mash') is not None and config['input_type'] != 'fasta' else [],
        INFORMS_nextclade = nextclade3.OUTPUT_INFORMS if 'nextclade' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        working_dir = config['working_dir']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.export_command_section(input, Path(output.HTML), params.working_dir)

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_reference_selection = refselection.OUTPUT_REPORT if config['fasta_ref'] is None and config['input_type'] != 'fasta' else refselection.OUTPUT_REPORT_EMPTY,
        report_preprocess_ampligone = preprocess.OUTPUT_AMPLIGONE_REPORT if 'ampligone' in config['analyses'] and config['input_type'] != 'fasta' else preprocess.OUTPUT_AMPLIGONE_REPORT_EMPTY,
        report_preprocess_clipping = preprocess.OUTPUT_CLIPPING_REPORT if 'ampligone' in config['analyses'] and config['input_type'] != 'fasta' else preprocess.OUTPUT_CLIPPING_REPORT_EMPTY,
        report_preprocess = preprocess.OUTPUT_REPORT if config['input_type'] != 'fasta' else [],
        report_iterative_mapping = iterativemapping.OUTPUT_REPORT if config['input_type'] != 'fasta' else [],
        report_nexclade_subtype = nextclade3.OUTPUT_SUBTYPE_REPORT if (config['nextclade'].get('db') is None) and ('nextclade' in config['analyses']) else nextclade3.OUTPUT_SUBTYPE_REPORT_EMPTY,
        report_nextclade = nextclade3.OUTPUT_REPORT if 'nextclade' in config['analyses'] else nextclade3.OUTPUT_REPORT_EMPTY,
        report_multi_allelic = multiallelicsites.OUTPUT_REPORT if config['input_type'] != 'fasta' else [],
        report_antivirals = antivirals.OUTPUT_REPORT if 'antivirals' in config['analyses'] else antivirals.OUTPUT_REPORT_EMPTY,
        report_commands = 'report/html-commands.iob',
        report_citations = core.OUTPUT_HTML_CITATIONS
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        ref_genome = config['fasta_ref'],
        input_type = config['input_type'],
        config_input = config['input'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        extra_info = [('Reference genome', Path(params.ref_genome).name if params.ref_genome is not None else '-')]
        input_file_str = ', '.join(f['name'] for _, input_files in params.config_input.items() for f in input_files)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            params.sample_name,
            datetime.datetime.now(),
            params.pipeline_info['version'], input_file_str, input_type=params.input_type, extra_data=extra_info))

        # Other sections
        report_structure = []

        # Core sections (shared)
        ReportPipeline.add_content_scrubbing(
            report_structure, params.input_type, input.reports_scrubbing)
        ReportPipeline.add_content_trim_basic_qc(
            report_structure, params.input_type, input.reports_downsampling, input.reports_trimming)
        ReportPipeline.add_content_contamination_check(
            report_structure, params.input_type, input.reports_contamination, None)

        # Add output sections
        if params.input_type != 'fasta':
            report_structure.extend([
                ('Reference selection', 'ref_selection', [Path(input.report_reference_selection)]),
                ('Pre-processing', 'pre_process', [Path(x) for x in (
                    input.report_preprocess_ampligone, input.report_preprocess_clipping, input.report_preprocess)]),
                ('Consensus extraction', 'consensus', [Path(input.report_iterative_mapping)]),
                ('Multi-allelic sites', 'multi_allelic', [Path(input.report_multi_allelic)]),
            ])
        report_structure.extend([
            ('Nextclade', 'nextclade', [Path(input.report_nexclade_subtype), Path(input.report_nextclade)]),
            ('Antiviral resistance', 'antiviral', [Path(input.report_antivirals)]),
            ('Commands', 'commands', [Path(input.report_commands)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
        ])
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    In this rule all summary files are combined into a complete summary output file.
    """
    input:
        core.OUTPUT_SUMMARY_INIT,
        lambda wildcards: human_read_scrubbing.get_summaries(config, wildcards.ext),
        lambda wildcards: downsampling.get_summaries(config, wildcards.ext),
        trimming.get_summaries(config),
        lambda wildcards: contamination_check_kraken.get_summaries(config, wildcards.ext),
        refselection.OUTPUT_SUMMARY if config['fasta_ref'] is None and config['input_type'] != 'fasta' else [],
        preprocess.OUTPUT_SUMMARY if config['input_type'] != 'fasta' else [],
        iterativemapping.OUTPUT_SUMMARY if config['input_type'] != 'fasta' else [],
        multiallelicsites.OUTPUT_SUMMARY if config['input_type'] != 'fasta' else [],
        nextclade3.OUTPUT_SUMMARY if 'nextclade' in config['analyses'] else [],
        antivirals.OUTPUT_SUMMARY if 'antivirals' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.combine_summary_data(input, Path(output.FILE), str(params.ext))
