from pathlib import Path

from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import human_read_scrubbing

#######################
# Included Snakefiles #
#######################
include: human_read_scrubbing.SNAKEFILE

#########
# Rules #
#########
rule all:
    input:
        HTML = config['output']['html'],
        OUT = human_read_scrubbing.get_output_io(config)


rule link_scrubbing_input:
    """
    Creates the FASTQ input for the human read scrubbing step.
    """
    output:
        FASTQ = human_read_scrubbing.INPUT_FASTQ if not 'fasta' in config['input'] else [],
        FASTA = human_read_scrubbing.INPUT_FASTA if 'fasta' in config['input'] else []
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        if 'fasta' in config['input']:
            snakemakeutils.dump_object([ToolIOFile(Path(config['input']['fasta'][0]['path']))], Path(output.FASTA))
        elif 'fastq_pe' in config['input']:
            snakemakeutils.dump_object([ToolIOFile(Path(f['path'])) for f in config['input']['fastq_pe']], Path(output.FASTQ))
        elif 'fastq_se' in config['input']:
            snakemakeutils.dump_object([ToolIOFile(Path(f['path'])) for f in config['input']['fastq_se']], Path(output.FASTQ))


rule report_command_section:
    """
    Creates the report section containing the tool commands.
    """
    input:
        INFORMS_scrubbing = human_read_scrubbing.get_command_informs(config)
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
        report_commands = rules.report_command_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        input_dict = config['input'],
        output_dir = config['output']['dir'],
        pipeline_info = config['script_info'],
        citation_keys = config['citations']
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
            key_citation=params.citation_keys['main'],
        ))

        # Add report content
        report_structure = []
        basepipeutils.add_content_scrubbing(
            report_structure, script_input.type_.value, input.reports_scrubbing)
        report_structure.extend([('Commands', 'commands', [Path(input.report_commands)])])
        snakepipelineutils.add_report_content(report, report_structure)
