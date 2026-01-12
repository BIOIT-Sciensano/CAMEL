import json
import re
from datetime import datetime
from importlib.resources import files
from pathlib import Path
from typing import Any, Optional, Union

import yaml

from camel.app.core.command import Command
from camel.app.core.reports.htmlcitation import HtmlCitation
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.config import config
from camel.app.core.errors import SnakemakeExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.snakemake import snakemakeutils
from camel.app.loggers import logger


def init_pipeline_report(output_path: Path, output_dir: Path, pipeline_info: dict[str, str]) -> HtmlReport:
    """
    Initializes an empty pipeline report.
    :param output_path: Output path
    :param output_dir: Output directory
    :param pipeline_info: Pipeline information
    :return: Report object
    """
    path_css = Path(str(files('camel').joinpath('resources/reports/style.css')))
    path_jquery = Path(str(files('camel').joinpath('resources/reports/jquery-3.2.1.min.js')))
    report = HtmlReport(output_path, output_dir, [path_jquery])
    report.initialize(pipeline_info['name'], path_css)
    report.add_pipeline_header(pipeline_info['title'])
    return report


def __get_failed_rule(stderr: str) -> Union[str, None]:
    """
    Returns the name of the rule that failed during Snakemake execution.
    :return: Name of the failed rule
    """
    for line in reversed(stderr.splitlines()):
        m = re.match(r'Error in rule (\w+):', line.strip())
        if not m:
            continue
        return m.group(1)
    return None


def run_snakemake(
        snakefile: str | Path, config_path: str | Path, targets: list[Path], working_dir: Path, threads: int = 8,
        resources: Optional[dict[str, Any]] = None, slurm_args: Optional[dict[str, int]] = None) -> Command:
    """
    Helper function to run snakemake workflows.
    :param snakefile: Workflow snakefile
    :param config_path: Path to configuration file
    :param targets: Target output files
    :param working_dir: Working directory
    :param threads: Number of threads to use
    :param resources: Dictionary of resources by keyword
    :param slurm_args: Dictionary of slurm arguments
    :return: None
    """
    if not working_dir.exists():
        working_dir.mkdir(parents=True)

    # Construct the base command
    command_parts = [
        'snakemake',
        *[str(x) for x in targets],
        '--snakefile', str(snakefile),
        '--configfile', str(config_path),
        '--cores', str(threads)
    ]

    # Add resources if they are specified
    if resources is not None:
        command_parts.append('--resources')
        for key, value in resources.items():
            command_parts.append(f'{key}={value}')

    # Add slurm submit file and parameters if specified
    if slurm_args is not None:
        command_parts.append(f'--cluster "{slurm_args["cluster"]}"')
        for key, value in slurm_args.items():
            if key != 'cluster':
                command_parts.append(f'--{key} {value}')

    # Create and run command
    command = Command(' '.join(command_parts))
    command.run(working_dir)
    if command.returncode != 0:
        rule_failed = __get_failed_rule(command.stderr)
        logger.error(f"Failed at rule: {rule_failed if rule_failed is not None else 'n/a'}")
        raise SnakemakeExecutionError(command.stdout, command.stderr, rule_failed)
    return command

def create_input_section(
        sample_name: str, date: datetime, pipeline_version: str, input_files: str, input_type: str,
        extra_data: Optional[list[tuple[str, str]]] = None, key_citation: str = None) -> HtmlReportSection:
    """
    Creates the input section for the HTML report.
    :param sample_name: Sample name
    :param date: Analysis date
    :param pipeline_version: Pipeline version
    :param input_files: Input files
    :param input_type: Input type
    :param extra_data: Extra data to include in the input section
    :param key_citation: Citation for the pipeline.
    :return: Input report section
    """
    table_data = [
        ['Sample:', sample_name],
        ['Analysis date:', date.strftime(config.date_fmt)],
        ['Pipeline version:', pipeline_version],
        ['Input files:', input_files],
        ['Input type:', input_type]
    ]
    if extra_data is not None:
        for key, value in extra_data:
            table_data.append([f'{key}:', value])
    section = HtmlReportSection('Input')
    section.add_table(table_data, table_attributes=[('class', 'information')])
    if key_citation is not None:
        section.add_header('Disclaimer', 2)
        section.add_paragraph('If you use this pipeline for your scientific work, please cite:')
        section.add_html_object(HtmlCitation.parse_from_json(key_citation))
    if input_type == 'fasta':
        section.add_warning_message(
            'The input file is in FASTA format, which is not compatible with some of the QC checks. '
            'Care should be taken when interpreting the results.')
    return section

def add_report_content(report: HtmlReport, report_structure: list[tuple[str, str, list[Path]]]) -> None:
    """
    Adds the content to the HTML report.
    :param report: Report
    :param report_structure: Report structure
    :return: None
    """
    # Add the overview section
    report.add_module_header('Sections')
    section = HtmlReportSection(None)
    overview_list = HtmlElement('ul')
    for title, key, _ in report_structure:
        list_item = HtmlElement('li')
        list_item.add_html_object(HtmlElement('a', title, [('href', f'#{key}')]))
        overview_list.add_html_object(list_item)
    section.add_html_object(overview_list)
    report.add_html_object(section)

    # Add the different sections
    for title, key, items in report_structure:
        report.add_module_header(title, key)
        for pickle in items:
            if not pickle.exists():
                continue
            section = snakemakeutils.load_object(pickle)[0].value
            report.add_html_object(section)
            section.copy_files(report.output_dir)
    report.save()

def combine_summary_data(snake_in: Any, path_out: Path, ext: str) -> None:
    """
    Combines the summary data into a single file.
    :param snake_in: Snakemake input
    :param path_out: Output path
    :param ext: Summary format (TSV / JSON)
    :return: None
    """
    if ext == 'json':
        json_dict = {}
        for summary_input in snake_in:
            with Path(summary_input).open() as handle_in:
                try:
                    json_dict = {**json_dict, **json.load(handle_in)}
                except json.JSONDecodeError as err:
                    logger.error(f"Could not parse {summary_input}")
                    raise err
        with path_out.open('w') as handle_out:
            json.dump(json_dict, handle_out, indent=2)
    elif ext == 'tsv':
        with path_out.open('w') as handle_out:
            for summary_input in snake_in:
                with open(summary_input) as handle_in:
                    handle_out.write(handle_in.read())
    else:
        raise ValueError(f"Invalid 'ext' value: {ext}")


def extract_fq_input(io_dict: Path, key_pe: Optional[str] = 'FASTQ_PE', key_se: Optional[str] = None,
                      keys_se: Optional[list[str]] = None, drop_empty: bool = False, read_type: str = 'PE') -> \
        dict[str, list[ToolIOFile]]:
    """
    Extracts a specific FASTQ input dictionary from the standardized FASTQ dictionary.
    :param io_dict: Path to input IO file
    :param key_pe: Key for paired end FASTQ files
    :param key_se: Key for single end FASTQ files
    :param keys_se: Separate keys for the forward and reverse SE FASTQ files
    :param drop_empty: If True, keys with no reads are dropped from the output
    :param read_type: Type of reads ('PE' or 'SE')
    :return: Reformatted dictionary
    """
    io = snakemakeutils.load_object(io_dict)
    output_dict = {}

    # Single end reads (no paired / orphaned reads available)
    if read_type == 'SE':
        output_dict[key_se] = io['SE']
        return output_dict

    # PE reads
    output_dict[key_pe] = io['PE']

    # Add SE reads
    if keys_se is not None:
        for key_orig, key_new in zip(['SE_FWD', 'SE_REV'], keys_se):
            try:
                output_dict[key_new] = io[key_orig]
            except KeyError:
                logger.warning(f"No '{key_orig}' input found")
    elif key_se is not None:
        se_reads = io.get('SE_FWD', []) + io.get('SE_REV', [])
        output_dict[key_se] = se_reads
    else:
        logger.debug("No key(s) provided for SE reads")

    # Remove keys that are empty
    if drop_empty:
        for key in list(output_dict.keys()):
            if (output_dict[key] is not None) and (len(output_dict[key]) > 0) and not (len(output_dict[key]) == 1 and output_dict[key][0].size == 0):
                continue
            logger.debug(f'Removing empty input: {key}')
            output_dict.pop(key)

    # Return the reformatted dictionary
    return output_dict

def create_empty_report_section(title: str, output_file: Path, header_level: int = 3) -> None:
    """
    Creates an empty report section.
    :param title: Section title
    :param output_file: Output file
    :param header_level: Header level
    :return: None
    """
    section = HtmlReportSection(title, header_level)
    section.add_paragraph('Analysis disabled.')
    snakemakeutils.dump_object([ToolIOValue(section)], output_file)

def generate_config_file(config_data: dict[str, Any], output_dir: Path, output_basename: str = 'config.yml') -> str:
    """
    Generates a configuration file for Snakemake in YAML file format.
    :param config_data: Configuration data
    :param output_dir: Output directory
    :param output_basename: Output basename
    :return: Path to config file
    """
    config_path = output_dir / output_basename
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle)
    logger.info(f"Configuration file created: {config_path}")
    return str(config_path)
