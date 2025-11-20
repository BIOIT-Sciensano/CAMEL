import json

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class CharacterizeNeisseriaCapsuleReporter(Tool):
    """
    Creates an HTML output report for the characterize_neisseria_capsule tool.
    """

    COLS_GENES = [
        {'key': 'allele_name', 'name': 'Gene', 'fmt': lambda x: f'<i>{x}</i>'},
        {'key': 'identity', 'name': '% identity', 'fmt': lambda x: f'{x:.2f}'},
        {'key': 'cov', 'name': '% covered', 'fmt': lambda x: f'{x*100:.2f}'},
        {'key': 'contig', 'name': 'Contig'},
        {'key': 'qstart', 'name': 'Start'},
        {'key': 'qend', 'name': 'End'}
    ]

    def __init__(self) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('CharacterizeNeisseriaCapsule reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("Serogrouping tool tabular output is required ('TSV')")
        if 'JSON' not in self._tool_inputs:
            raise InvalidToolInputError("Serogrouping tool JSON output is required ('JSON')")
        if 'detector' not in self._input_informs:
            raise InvalidToolInputError("Detector informs are required")
        super()._check_input()

    def __add_overview_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the detected serogroup to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Rename columns
        header = {'SG': 'Predicted serogroup', 'Notes': 'Notes'}

        # Create table data
        table_data = [[r['SG'], r['Notes']] for r in data.to_dict('records')]
        section.add_table(table_data, list(header.values()), [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection(
            'Capsule characterization', subtitle=self._input_informs['detector']['_name_full'])

        # Add overview table
        section.add_header('Overview', 3)
        self.__add_overview_table(section)

        # Parse detected genes
        section.add_header('Detected genes', 3)
        with open(self._tool_inputs['JSON'][0].path) as handle:
            json_info = json.load(handle)
        data_genes = pd.DataFrame(json_info['Serogroup'][0]['genes'])

        # Add the detected genes to the report (if there are any)
        if len(data_genes) == 0:
            section.add_paragraph('No capsule genes detected.')
        else:
            data_genes['color'] = data_genes.apply(
                lambda x: CharacterizeNeisseriaCapsuleReporter.__get_row_color(x), axis=1)
            data_genes.sort_values(by='allele_name', inplace=True)

            # Add table with overview of the detected genes
            table_data = [
                [HtmlTableCell(col.get('fmt', str)(row[col['key']]), color=row['color']) for col in
                 CharacterizeNeisseriaCapsuleReporter.COLS_GENES] for row in data_genes.to_dict('records')
            ]
            section.add_table(
                table_data, [c['name'] for c in CharacterizeNeisseriaCapsuleReporter.COLS_GENES], [('class', 'data')])

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def __get_row_color(row: pd.Series) -> str:
        """
        Returns the row color for a detected gene.
        :param row: Input row
        :return: Color
        """
        if row['cov'] == 1.0 and row['identity'] == 100.0:
            return 'green'
        elif row['cov'] == 1.0:
            return 'lightgreen'
        return 'grey'
