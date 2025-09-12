from pathlib import Path
from typing import Any

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class BacMetReporter(Tool):
    """
    Creates an HTML report for the BacMet analysis.
    """

    COLUMN_MAPPING = {
        'BacMet_ID': {'title': 'ID'},
        'Gene_name': {'title': 'Gene', 'fmt': lambda x: f'<i>{x}</i>'},
        'slen': {'title': 'Length (AA)', 'fmt': lambda x: f'{x:,}'},
        'pident': {'title': 'Identity (%)', 'fmt': lambda x: f'{x:.2f}'},
        'perc_covered': {'title': 'Covered (%)', 'fmt': lambda x: f'{x:.2f}'},
        'Accession': {'title': 'Accession'},
        'Organism': {'title': 'Organism', 'fmt': lambda x: f'<i>{x}</i>'},
        'Compound': {'title': 'Compound'}
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('BacMet reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('TSV input is required')
        if 'blastp' not in self._input_informs:
            raise InvalidToolInputError("blastp informs input is required")
        if 'filtering' not in self._input_informs:
            raise InvalidToolInputError("filtering informs input is required")
        super()._check_input()

    @staticmethod
    def __get_row_color(row: dict[str, Any]) -> str:
        """
        Returns the color for the given row.
        :param row: Input row
        :return: Row color
        """
        if row['pident'] == 100.0 and row['perc_covered'] == 100.0:
            return 'green'
        if row['pident'] < 100.0 and row['perc_covered'] == 100.0:
            return 'lightgreen'
        return 'grey'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('BacMet', subtitle=self._input_informs['blastp']['_name'])

        # Parameters
        section.add_table([
            ['Min % identity (AA):', f"{self._input_informs['filtering']['params']['min_id']:.2f}%"],
            ['Min % covered (AA):', f"{self._input_informs['filtering']['params']['min_cov']:.2f}%"],
            ['Database:', 'Experimentally confirmed genes']
        ], None, [('class', 'information')])

        # Results table
        data_in = pd.read_table(self._tool_inputs['TSV'][0].path)
        table_data = [
            [HtmlTableCell(d.get('fmt', lambda x: x)(row[col]), color=BacMetReporter.__get_row_color(row)) for
             col, d in BacMetReporter.COLUMN_MAPPING.items()] for row in data_in.to_dict('records')
        ]
        section.add_table(table_data, [c['title'] for c in BacMetReporter.COLUMN_MAPPING.values()], [('class', 'data')])

        # Download link
        relative_path = Path('bacmet', 'blastp.tsv')
        section.add_link_to_file('Download (TSV)', relative_path)
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path)

        # Last update data
        path_version = self._tool_inputs['DB'][0].path / 'VERSION'
        with path_version.open() as handle:
            version = handle.readline().strip()
        section.add_paragraph(f'Database version: {version}')

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
