from pathlib import Path

from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.snakefiles import abritamr, assembly, core

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
        HTML = config['output']['html'],
        TSV = config['output']['tsv']

rule report_pickle_citations:
    """
    This rule creates a pickle with a report section containing the citations.
    """
    output:
        HTML = Path(config['working_dir']) / 'report' / 'html-citations.iob'
    params:
        citation_keys = config['citations']
    run:
        from camel.app.core.reports import reportutils
        section = reportutils.create_citations_section(
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
        from camel.app.scriptutils.basepipe import basepipeutils
        basepipeutils.export_command_section(input, Path(output.HTML), Path(params.dir_))

rule report_combine_all:
    """
    Rule to combine report sections into a single output report.
    """
    input:
        report_abritamr = abritamr.OUTPUT_REPORT,
        report_citations = rules.report_pickle_citations.output.HTML,
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        sample_name = config['input']['sample_name'],
        input_dict = config['input'],
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        citation_keys = config['citations']
    run:
        import datetime
        from camel.app.scriptutils import model

        # Add header section
        report = snakepipelineutils.init_pipeline_report(
            Path(output.HTML), Path(params.output_dir), params.pipeline_info)
        report.add_html_object(snakepipelineutils.create_input_section(
            sample_name=params.input_dict['sample_name'],
            date=datetime.datetime.now(),
            pipeline_version=params.pipeline_info['version'],
            input_files=params.input_dict['fasta']['name'],
            input_type=model.InputType.FASTA.value,
        ))

        # Add report content
        report_structure = [
            ('AbriTAMR', 'abrit', [Path(input.report_abritamr)]),
            ('Citations', 'citations', [Path(input.report_citations)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        snakepipelineutils.add_report_content(report, report_structure)

rule summary_combine_all:
    """
    This rule copies the resources output summary to the config output summary.
    """
    input:
        TSV = str(abritamr.OUTPUT_SUMMARY).format(ext='tsv')
    output:
        TSV = config['output']['tsv']
    shell:
        """
        cp {input.TSV} {output.TSV}
        """
