import pandas as pd
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class MenDeVARReporter(Tool):
    """
    Creates an HTML output report for the MenDeVAR indexes.
    """

    def __init__(self) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('MenDeVAR reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'mendevar' not in self._input_informs:
            raise InvalidToolInputError("MenDeVAR informs are required ('mendevar')")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("MenDeVAR tabular output is required ('TSV')")
        super()._check_input()

    def __add_bexsero_index_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the detected alleles and their respective MenDeVAR indexes to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data_hits = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove columns with Trumenba indexes
        bexsero_data = data_hits.loc[:, ~data_hits.columns.str.startswith('Trumenba')]

        # Create table data
        table_data = []
        for values in bexsero_data.fillna('-').itertuples(index=False, name=None):
            row = [v if not isinstance(v, float) else f'{v:.2f}' for v in values]

            # Color the cell with the allele status
            bexsero_status = list(bexsero_data.columns).index('Bexsero allele status')
            row[bexsero_status] = HtmlTableCell(row[bexsero_status], MenDeVARReporter.__get_color(row[bexsero_status]))
            table_data.append(row)

        section.add_table(table_data, list(bexsero_data.columns), [('class', 'data')])
        section.add_paragraph(
            '<b>Note:</b> Imperfect matches are marked with an asterisk (*) and are not used for determining the '
            'MenDeVAR index.')

    def __add_trumenba_index_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the detected alleles and their respective MenDeVAR indexes to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data_hits = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove columns with Bexsero indexes
        trumenba_data = data_hits.loc[:, ~data_hits.columns.str.startswith('Bexsero')]

        # Create table data
        table_data = []
        for values in trumenba_data.fillna('-').itertuples(index=False, name=None):
            row = [v if not isinstance(v, float) else f'{v:.2f}' for v in values]

            # Color the cell with the allele status
            trumenba_status = list(trumenba_data.columns).index('Trumenba allele status')
            row[trumenba_status] = HtmlTableCell(
                row[trumenba_status], MenDeVARReporter.__get_color(row[trumenba_status]))
            table_data.append(row)

        section.add_table(table_data, list(trumenba_data.columns), [('class', 'data')])
        section.add_paragraph(
            '<b>Note:</b> Imperfect matches are marked with an asterisk (*) and are not used for determining the '
            'MenDeVAR index.')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('MenDeVAR indexes')

        # Create overview table for the MenDeVAR Bexsero indexes
        section.add_header('Bexsero MenDeVAR index', 3)
        bexsero_status = self._input_informs['mendevar']['bexsero_status']
        section.add_table(
            [[HtmlTableCell(bexsero_status, MenDeVARReporter.__get_color(bexsero_status))]], ['Bexsero status'],
            [('class', 'data')])

        # Create table with Bexsero indexes
        self.__add_bexsero_index_table(section)
        section.add_horizontal_line()

        # Create overview table for the MenDeVAR Trumenba indexes
        section.add_header('Trumenba MenDeVAR index', 3)
        trumenba_status = self._input_informs['mendevar']['trumenba_status']
        section.add_table(
            [[HtmlTableCell(trumenba_status, MenDeVARReporter.__get_color(trumenba_status))]], ['Trumenba status'],
            [('class', 'data')])

        # Create table with Trumenba indexes
        self.__add_trumenba_index_table(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def __get_color(status: str) -> str:
        """
        Returns the color for the corresponding MenDeVAR index.
        :param status: MenDeVAR index
        :return: Status cell
        """
        if status == 'exact match':
            return 'green'
        if status == 'cross-reactive':
            return 'yellow'
        if status == 'none':
            return 'red'
        return 'grey'
