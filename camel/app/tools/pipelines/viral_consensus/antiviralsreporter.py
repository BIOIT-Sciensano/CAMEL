import json

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.reports.htmltableformatter import HtmlTableFormatter, FormatEntry
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


def get_color(inhibition: str) -> HtmlTableCell:
    """
    Returns a colored cell for the target inhibition level.
    :return: Formated table cell
    """
    if inhibition == 'NI':
        return HtmlTableCell(inhibition, color='green')
    return HtmlTableCell(inhibition, color='red')


class AntiviralsReporter(Tool):
    """
    Reporter class for the antiviral mutation detection.
    """
    COLS_MUTS: list[FormatEntry] = [
        {'key': 'subtype', 'title': 'Subtype'},
        {'key': 'segment', 'title': 'Segment'},
        {'key': 'type', 'title': 'Category'},
        {'key': 'mutation', 'title': 'Mutation'},
    ]

    COLS_ASSOCIATIONS: list[FormatEntry] = [
        {'key': 'category', 'title': 'Category'},
        {'key': 'key', 'title': 'Key'},
        {'key': 'antiviral', 'title': 'Antiviral'},
        {'key': 'resistance', 'title': 'Resistance', 'fmt': get_color},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('antiviral reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError('JSON input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        with self._tool_inputs['JSON'][0].path.open() as handle:
            data_detection = json.load(handle)

        # Initialize the report section
        section = HtmlReportSection('Antiviral resistance')

        # Add table with mutations
        section.add_header('Detected mutations', 4)
        if len(data_detection['mutations']) == 0:
            section.add_paragraph('No mutations associated with antiviral resistance detected.')
        else:
            section.add_table(
                HtmlTableFormatter.format_table_data(
                    pd.DataFrame(data_detection['mutations']), AntiviralsReporter.COLS_MUTS),
                [col['title'] for col in AntiviralsReporter.COLS_MUTS], [('class', 'data')]
            )
            section.add_paragraph(
                'Note: Only mutations from the database, associated with antivirals, are shown here. A complete list '
                'of all mutations is included in the Nextclade output.')

        # Add table with associations
        section.add_header('Associations', 4)
        if len(data_detection['associations']) == 0:
            section.add_paragraph('No antiviral resistance associations detected.')
        else:
            section.add_table(
                HtmlTableFormatter.format_table_data(
                    pd.DataFrame(data_detection['associations']), AntiviralsReporter.COLS_ASSOCIATIONS),
                [col['title'] for col in AntiviralsReporter.COLS_ASSOCIATIONS], [('class', 'data')]
            )

        # Add a paragraph with additional information
        section.add_header('Additional information', 4)
        section.add_paragraph(
            "The mutations are extracted from the Nextclade output. Abbreviations: normal inhibition (<b>NI</b>), "
            "reduced inhibition (<b>RI</b>), highly reduced inhibition (<b>HRI</b>).")

        # Save output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
