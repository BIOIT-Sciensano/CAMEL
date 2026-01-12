from datetime import datetime
from importlib.resources import files
from pathlib import Path
from collections.abc import Sequence
from typing import Any

from camel.app.config import config
from camel.app.core.reports.htmlcitation import HtmlCitation
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.loggers import logger


def init_report(path_out: Path, key: str, title: str, dir_out: Path | None = None) -> HtmlReport:
    """
    Initializes an empty HTML report.
    :param path_out: Report output path
    :param key: Report key
    :param title: Report title (can contain additional formatting)
    :param dir_out: Output directory (defaults to parent directory of the report)
    """
    # Resources
    path_css = Path(str(files('camel').joinpath('resources/reports/style.css')))
    path_jquery = Path(str(files('camel').joinpath('resources/reports/jquery-3.2.1.min.js')))

    # Report object
    if not path_out.parent.is_dir():
        logger.info(f"Creating directory to store report: {path_out}")
        path_out.parent.mkdir(parents=True)
    report = HtmlReport(path_out, dir_out if dir_out else path_out.parent, [path_jquery])
    report.output_dir.mkdir(exist_ok=True, parents=True)
    report.initialize(key, path_css)
    report.add_pipeline_header(title)
    return report


def create_overview_section(
        version: str,
        dataset_name: str,
        input_file_str: str,
        date: datetime | None = None,
        input_type: str | None = None,
        extra_data: list[Sequence[str]] | None = None,
        key_citation: str | None = None,
        warnings: list[str] | None = None
    ) -> HtmlReportSection:
    """
    Creates the overview section for the HTML output report.
    :param dataset_name: Dataset name
    :param input_file_str: String of the input files
    :param date: Analysis date
    :param version: Pipeline or tool version
    :param input_type: Input type
    :param extra_data: Extra data to include in the input section
    :param key_citation: Citation for the pipeline.
    :param warnings: List of warning to display in the report
    :return: Input report section
    """
    # Main information
    table_data = [
        ['Dataset:', dataset_name],
        ['Analysis date:', (date if not date is None else datetime.now()).strftime(config.date_fmt)],
        ['Version:', version],
        ['Input files:', input_file_str],
    ]

    # Additional info
    if input_type is not None:
        table_data.append(['Input type:', input_type])
    if extra_data is not None:
        for key, value in extra_data:
            table_data.append([f'{key}:', value])

    # Create and return the section
    section = HtmlReportSection('Input')
    section.add_table(table_data, table_attributes=[('class', 'information')])

    # Add the citation (optional)
    if key_citation is not None:
        section.add_header('Disclaimer', 2)
        section.add_paragraph('If you use this for your scientific work, please cite:')
        section.add_html_object(HtmlCitation.parse_from_json(key_citation))

    # Add warning messages (optional)
    if warnings is not None:
        for warning_str in warnings:
            section.add_warning_message(warning_str)
    return section

def create_commands_section(
        tool_informs: list[dict[str, Any]],
        dir_: Path) -> HtmlReportSection:
    """
    Creates a section with an overview of the commands.
    :param tool_informs: Tool informs
    :param dir_: Working directory
    :return: Commands section
    """
    section = HtmlReportSection('Commands')
    logger.debug(f"Exporting command for {len(tool_informs)} tools")
    for informs in tool_informs:
        header = f"{informs['_name_full']} - {informs['_tag']}" if '_tag' in informs else informs['_name']
        section.add_header(header, 3)
        command_txt = informs['_command'].replace(str(dir_), '$WORKING')
        command_txt = command_txt.replace('\n', '<br />\n')
        section.add_html_object(HtmlElement('code', command_txt, [('class', 'command')]))
    return section


def create_citations_section(keys_other: list[str], key_main: str | None = None) -> HtmlReportSection:
    """
    Creates the report section with the citations.
    :param keys_other: List of key for citations for tools and databases
    :param key_main: Key for the main citation of the workflow
    :return: Citations section
    """
    section_citations = HtmlReportSection('Citations')
    if key_main is not None:
        section_citations.add_header('Pipeline', 3)
        section_citations.add_html_object(HtmlCitation.parse_from_json(key_main))
    section_citations.add_header('Tools and databases', 3)
    for citation_key in keys_other:
        section_citations.add_html_object(HtmlCitation.parse_from_json(citation_key))
    return section_citations
