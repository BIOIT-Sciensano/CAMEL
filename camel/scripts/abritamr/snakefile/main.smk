from pathlib import Path

from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import abritamr, assembly, core

#######################
# Included Snakefiles #
#######################
include: core.SNAKEFILE
include: abritamr.SNAKEFILE
include: assembly.SNAKEFILE

#########
# Rules #
#########
rule all:
    input:
        HTML = config['output_report'],
        TSV = config['output_tabular']


rule report_pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-citations.io'
    params:
        citation_keys = config['citations']
    run:
        section = SnakePipelineUtils.create_citations_section(
            params.citation_keys['other'], params.citation_keys['main'])
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule report_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_abritamr_run =  abritamr.OUTPUT_RUN_INFORMS,
        INFORMS_abritamr_report = abritamr.OUTPUT_REPORT_REPORT_INFORMS
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
        report_abritamr = abritamr.OUTPUT_REPORT,
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output_report']
    params:
        sample_name = config['sample_name'],
        input_dict = config['input'],
        input_type = config['input_type'],
        output_dir = config['output_dir'],
        pipeline_info = config['pipeline'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline

        # Add header section
        report = SnakePipelineUtils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(SnakePipelineUtils.create_input_section(
            sample_name=params.sample_name,
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=ReportPipeline.format_input_string(params.input_dict),
            input_type=params.input_type,
            key_citation=params.citation_keys['main']
        ))

        # Add report content
        report_structure = [
            ('AbriTAMR', 'abrit', [Path(input.report_abritamr)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        SnakePipelineUtils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    This rule copies the resources output summary to the config output summary.
    """
    input:
        TSV = abritamr.OUTPUT_SUMMARY
    output:
        TSV = config['output_tabular']
    shell:
        """
        cp {input.TSV} {output.TSV}
        """
