import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.reports.htmltableformatter import HtmlTableFormatter, FormatEntry
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ShigaTyperReporter(Tool):
    """
    Creates an HTML output report for the ShigaTyper tool.
    """

    COLS: list[FormatEntry] = [
        {'key': 'Hit', 'title': 'Hit'},
        {'key': 'Number of reads', 'title': 'Nb. of reads', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'Length Covered', 'title': 'Length covered', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'reference length', 'title': 'Reference length', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': '% covered', 'title': '% covered', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'Number of variants', 'title': 'Nb. of variants', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': '% accuracy', 'title': '% accuracy', 'fmt': HtmlTableFormatter.FLOAT_FMT}
    ]

    def __init__(self) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('ShigaTyper reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("ShigaTyper tabular output is required ('TSV')")
        if 'shigatyper' not in self._input_informs:
            raise InvalidToolInputError("ShigaTyper informs are required")
        super()._check_input()

    def __add_main_output(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the ShigaTyper output to the report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input files
        main_output = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove the sample ID column
        main_output.pop('sample')
        main_output.fillna('n/a', inplace=True)

        # Create table data
        main_table = []
        for values in main_output.itertuples(index=False, name=None):
            row = list(values)
            main_table.append(row)

        # Rename columns
        header = ['Prediction',	'ipaB', 'Notes']

        section.add_table(main_table, header, [('class', 'data')])

    def __add_shigatyper_hits(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the ShigaTyper hits details to the report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input files
        gene_hits = pd.read_table(self._tool_inputs['TSV_HITS'][0].path)
        gene_hits['color'] = gene_hits.apply(lambda x: ShigaTyperReporter.get_color(x), axis=1)

        # Remove the sample ID column
        gene_hits.pop('Unnamed: 0')

        # Add potentially missing columns
        col_list = ['Hit', 'Number of reads', 'Length Covered', 'reference length',
                    '% covered', 'Number of variants', '% accuracy']
        missing_columns = [col for col in col_list if col not in gene_hits.columns]
        gene_hits[missing_columns] = '-'

        # Rename columns
        header = [c['title'] for c in ShigaTyperReporter.COLS]
        section.add_table(HtmlTableFormatter.format_table_data(
            gene_hits, ShigaTyperReporter.COLS, use_colors=True), header, [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection('ShigaTyper', subtitle=self._input_informs['shigatyper']['_name'])
        species = self._input_informs['shigatyper']['species']
        section.add_table(
            [[HtmlTableCell(species)]], ['Serotyping'], [('class', 'data')])

        # Create table with hits
        section.add_header('Overview', 3)
        self.__add_main_output(section)
        if species != 'Not Shigella or EIEC':
            self.__add_shigatyper_hits(section)

            # Add warning when pseudoreads were used for the execution of the tool
        if 'pseudo_reads' in self._parameters:
            section.add_warning_message("The tool is executed on simulated reads.")

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def get_color(row_in: pd.Series) -> str:
        """
        Colors the rows based on the statistics.
        :param row_in: Input row
        :return: Color of the row
        """
        if row_in['% covered'] == 100.0 and row_in['% accuracy'] == 100.0:
            return 'green'
        if row_in['% covered'] == 100.0:
            return 'lgreen'
        return 'grey'
