import json
from collections import defaultdict

import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlexpandablediv import HtmlExpandableDiv
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell
from camelcore.app.reports.htmltableformatter import FormatEntry, HtmlTableFormatter

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.tools.pipelines.viral_consensus.antiviralsdetection import AntiviralsDetection

# Resistance codes (influenza & RSV): code -> (color, description)
RESISTANCE_CODES: dict[str, tuple[str, str]] = {
    'NI': ('green', 'Normal inhibition, or described polymorphism without impact on sensitivity'),
    'NI/RI': ('yellow', 'Normal to reduced inhibition'),
    'RI': ('yellow', 'Reduced inhibition'),
    'NI/RI/HRI': ('yellow', 'Normal to highly reduced inhibition'),
    'PR': ('yellow', 'Mutation associated with possible resistance'),
    'RI/HRI': ('red', 'Reduced to highly reduced inhibition'),
    'HRI': ('red', 'Highly reduced inhibition'),
    'R': ('red', 'Mutation associated with resistance'),
}


def get_color(inhibition: str) -> HtmlTableCell:
    """
    Returns a colored cell for the target inhibition level.
    :return: Formated table cell
    """
    color, _ = RESISTANCE_CODES.get(inhibition, ('red', None))
    return HtmlTableCell(inhibition, color=color)


class AntiviralsReporter(Tool):
    """
    Reporter class for the antiviral mutation detection.
    """
    COLS_MUTS: list[FormatEntry] = [
        {'key': 'species', 'title': 'Species'},
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
        super().__init__('antiviral reporter', '0.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError('JSON input is required')
        super()._check_input()

    @staticmethod
    def _get_positions_not_assessable(mutations: list[dict]) -> list[list[str]]:
        """
        Groups the database mutations that could not be assessed by segment.
        :param mutations: Mutations that could not be assessed
        :return: Table rows with the segment and the positions
        """
        positions_by_segment = defaultdict(set)
        for mutation in mutations:
            position = AntiviralsDetection.parse_position(mutation['mutation'])
            if position is not None:
                positions_by_segment[mutation['segment']].add(position)
        return [
            [segment, ', '.join(str(p) for p in sorted(positions))]
            for segment, positions in sorted(positions_by_segment.items())
        ]

    def _execute_tool(self) -> None:
        """
        Executes the tool.
        :return: None
        """
        with self._tool_inputs['JSON'][0].path.open() as handle:
            data_detection = json.load(handle)

        # Initialize the report section
        section = HtmlReportSection('Antiviral resistance')
        muts_not_assessable = data_detection.get('mutations_not_assessable', [])
        if len(muts_not_assessable) > 0:
            section.add_alert(
                f'{len(muts_not_assessable)} mutation(s) from the database could not be assessed due to missing '
                f'sequence data (e.g. low coverage). Resistance cannot be excluded for these positions.', 'warning')

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

        # Add table with positions that could not be assessed
        if len(muts_not_assessable) > 0:
            section.add_header('Positions not assessed', 4)
            section.add_table(
                AntiviralsReporter._get_positions_not_assessable(muts_not_assessable),
                ['Segment', 'Positions'], [('class', 'data')])

        # Add a paragraph with additional information
        section.add_header('Additional information', 4)
        section.add_paragraph('The mutations are extracted from the Nextclade output.')
        div = HtmlExpandableDiv('antivirals-resistance-codes', 'resistance codes')
        div.add_table(
            data=[[get_color(code), description] for code, (_, description) in RESISTANCE_CODES.items()],
            column_names=['Resistance code', 'Description'],
            table_attributes=[('class', 'data')],
        )
        section.add_html_object(div)

        # Save output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
