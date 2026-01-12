import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class ShigEiFinderReporter(Tool):
    """
    Creates an HTML output report for the ShigEiFinder tool.
    """

    def __init__(self) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('ShigEiFinder reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("ShigEifinder tabular output is required ('TSV')")
        if 'shigeifinder' not in self._input_informs:
            raise InvalidToolInputError("ShigEiFinder informs are required")
        super()._check_input()

    def __add_shigeifinder_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the ShigEiFinder output to the report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data_hits = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove the sample ID column
        data_hits.pop('#SAMPLE')

        # Create table data
        table_data = []
        for values in data_hits.fillna('n/a').itertuples(index=False, name=None):
            row = list(values)
            table_data.append(row)

        # Rename columns
        header = ['ipaH', 'Virulence plasmid', 'Cluster', 'Serotype', 'O antigen', 'H antigen', 'Notes']

        section.add_table(table_data, header, [('class', 'data')])

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection('ShigEiFinder', subtitle=self._input_informs['shigeifinder']['_name_full'])
        species = self._input_informs['shigeifinder']['species']
        section.add_table(
            [[HtmlTableCell(f'<i>{species}</i>')]], ['Serotyping'], [('class', 'data')])

        # Create table with hits
        section.add_header('Overview', 3)
        self.__add_shigeifinder_table(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
