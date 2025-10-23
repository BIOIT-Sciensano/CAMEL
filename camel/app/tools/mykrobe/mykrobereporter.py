from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class MykrobeReporter(Tool):
    """
    Parses Mykrobe csv output reports and generates an HTML report section for the final report.
    """

    TITLE = 'Mykrobe'

    COLS = [
        {'key': 'level'},
        {'key': 'id'},
        {'key': 'percent', 'fmt': lambda x: f'{x:.2f}' if isinstance(x, float) else x},
        {'key': 'depth', 'fmt': lambda x: f'{int(x):,}' if isinstance(x, float) else x}
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mykrobe Reporter', '0.1')
        self._section = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'CSV' not in self._tool_inputs:
            raise InvalidToolInputError("Mykrobe output is required ('CSV')")
        if 'mykrobe' not in self._input_informs:
            raise InvalidToolInputError("Mykrobe informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection(MykrobeReporter.TITLE,
                                    subtitle=self._input_informs['mykrobe']['_name'])
        # Antibiotic sensitivity
        if 'show_amr' in self._parameters:
            section.add_header('Antibiotic susceptibility', 3)
            self.__add_antibiotic_sensitivity(section)
            section.add_html_object(
                HtmlElement('a', 'Description of the fields of the tables',
                            [('href', 'https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output')]))

        # Genotyping
        section.add_header(self._parameters['custom_header'].value, 3)
        group = self._input_informs['mykrobe']['phylo_group']
        section.add_table(
            [[HtmlTableCell(group)]], ['High level identification'], [('class', 'data')])
        self.__add_lineage_information(section)

        # Additional information
        self.__add_output_table_link(section)
        section.add_html_object(
            HtmlElement('a', 'Description of the fields of the tables',
                        [('href', 'https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output')]))
        self.__add_database_information(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __add_antibiotic_sensitivity(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the antibiotic sensitivity.
        :return: None
        """
        header = ['Antibiotic', 'Susceptibility', 'Variants', 'Genes']
        data = []
        results = pd.read_csv(self._tool_inputs['CSV'][0].path)

        # Replace all nan by dashes
        results.replace(np.nan, '-', inplace=True)
        for i in range(results.shape[0]):
            data.append([results.iloc[i, 1],
                         results.iloc[i, 2],
                         ', '.join(results.iloc[i, 3].split(';')),
                         ', '.join(results.iloc[i, 4].split(';'))])

        # Start writing in the report the table and the headers
        section.add_table(data, header, [('class', 'data')])

    def __add_lineage_information(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the lineage information.
        :return: None
        """
        results = pd.read_csv(self._tool_inputs['CSV'][0].path)
        results.replace(np.nan, '-', inplace=True)

        # Create table data
        levels = {'Phylogenetic group': 'phylo', 'Species': 'species', 'Lineage': 'lineage'}
        table_list = []

        for key, value in levels.items():
            output = {'level': key}
            row = results.loc[:, results.columns.str.startswith(value)].drop_duplicates()
            row.columns = ['id', 'percent', 'depth']
            output.update(row.iloc[0].to_dict())
            table_list.append([MykrobeReporter.__format_cell(
                output[col['key']], col) for col in MykrobeReporter.COLS])

        # Rename columns
        header = ['Taxonomic level', 'Identification', 'Percent of probe that\nhas any coverage',
                  'Average depth\nacross probe']
        section.add_table(table_list, header, [('class', 'data')])

    def __add_output_table_link(self, section: HtmlReportSection) -> None:
        """
        Adds link to the output table (tsv) for this assay.
        :return: None
        """
        relative_path = Path('mykrobe', 'summary_out.tsv')
        section.add_file(self._tool_inputs['CSV'][0].path, relative_path)
        section.add_link_to_file("Download (CSV)", relative_path)

    def __add_database_information(self, section: HtmlReportSection) -> None:
        """
        Adds the database information to the report.
        :return: None
        """
        section.add_paragraph('Last updated: {}'.format(self._input_informs['mykrobe'].get(
            'last_update_date', '{LAST_UPDATE_DATE}')))

    @staticmethod
    def __format_cell(value: Any, col: dict[str, Any]) -> HtmlTableCell:
        """
        Formats the corresponding table cell.
        :param value: Input value
        :param col: Column metadata
        :return: HTML table cell
        """
        if 'fmt' not in col:
            return HtmlTableCell(str(value))
        return HtmlTableCell(col['fmt'](value))
