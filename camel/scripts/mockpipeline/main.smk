from pathlib import Path

from camel.snakefiles import trimming_illumina, downsampling, trimming_ont, trimming, quast, \
    contamination_check_kraken, quality_checks, confindr, gene_detection, assembly, core, human_read_scrubbing, \
    read_simulation, variant_calling, variant_filtering
from camel.scripts.mycobacteriumpipeline.snakefile import snpit
from camel.app.core.snakemake import snakepipelineutils

#######################
# Included snakefiles #
#######################
include: core.SNAKEFILE
include: human_read_scrubbing.SNAKEFILE
include: read_simulation.SNAKEFILE
include: downsampling.SNAKEFILE
include: trimming_illumina.SNAKEFILE
include: trimming_ont.SNAKEFILE
include: assembly.SNAKEFILE
include: variant_calling.SNAKEFILE
include: variant_filtering.SNAKEFILE
include: quast.SNAKEFILE
include: contamination_check_kraken.SNAKEFILE
include: confindr.SNAKEFILE
include: quality_checks.SNAKEFILE
include: gene_detection.SNAKEFILE
include: snpit.SNAKEFILE

#########s
# Rules #
#########
rule all:
    """
    Ensures that the required output files are generated.
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
        INFORMS_simulation =  read_simulation.OUTPUT_INFORMS if config['input']['type'] == 'fasta' else [],
        INFORMS_downsampling = downsampling.get_command_informs(config),
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_assembly = assembly.get_command_informs(config),
        INFORMS_quast = quast.OUTPUT_INFORMS,
        INFORMS_busco = quast.OUTPUT_INFORMS_BUSCO,
        INFORMS_contamination = contamination_check_kraken.get_command_informs(config),
        INFORMS_confindr = confindr.get_command_informs(config),
        INFORMS_ncbi_amr = gene_detection.OUTPUT_INFORMS.format(db='ncbi_amr') if 'ncbi_amr' in config['analyses'] else [],
        INFORMS_snpit = snpit.OUTPUT_INFORMS if 'snpit' in config['analyses'] else []
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

rule report_create:
    """
    Creates the output HTML report.
    """
    input:
        reports_scrubbing = human_read_scrubbing.get_reports(config),
        reports_downsampling = downsampling.get_reports(config),
        reports_trimming = trimming.get_reports(config),
        report_quast = quast.OUTPUT_REPORT,
        reports_contamination = contamination_check_kraken.get_reports(config),
        report_confindr = confindr.get_report(config),
        report_adv_qc = quality_checks.OUTPUT_REPORT.format(input_type=config['input']['type']),
        report_ncbi_amr = gene_detection.get_gene_detection_report('ncbi_amr', config),
        report_snpit = snpit.OUTPUT_REPORT if 'snpit' in config['analyses'] else snpit.OUTPUT_REPORT_EMPTY,
        report_commands = rules.report_create_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        input_dict = config['input'],
        gene_detection_method = config['gene_detection']['options']['method']
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
            extra_data=[('Gene detection method', params.gene_detection_method)]
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
        report_structure.append(('Gene detection', 'gene_detection', [Path(input.report_ncbi_amr)]))
        report_structure.append(('Species identification', 'identification', [Path(input.report_snpit)])),
        report_structure.append(('Commands', 'commands', [Path(input.report_commands)]))

        # Export the report
        snakepipelineutils.add_report_content(report, report_structure)

rule summary_combine:
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
        snpit.OUTPUT_SUMMARY if 'snpit' in config['analyses'] else []
    output:
        FILE = 'summary/output.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.combine_summary_data(input, Path(output.FILE), str(params.ext))
