import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class GMatsReporter(Tool):
    """
    Creates an HTML output report for the gMATS assay.
    """

    def __init__(self) -> None:
        """
        Initializes the reporter.
        """
        super().__init__('gMATS reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'gmats' not in self._input_informs:
            raise InvalidToolInputError("gMATS informs are required ('gmats')")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("gMATS tabular output is required ('TSV')")
        super()._check_input()

    def __add_alleles_table(self, section: HtmlReportSection) -> None:
        """
        Adds a table with the detected alleles to the output report section.
        :param section: Report output section
        :return: None
        """
        # Parse the input file
        data_hits = pd.read_table(self._tool_inputs['TSV'][0].path)

        # Remove the last gMATS status column (already shown above)
        data_hits.pop('gMATS status')

        # Create table data
        table_data = []
        for values in data_hits.itertuples(index=False, name=None):
            row = list(values)

            # Color the cell with the allele status
            idx_status = list(data_hits.columns).index('Allele status')
            row[idx_status] = HtmlTableCell(row[idx_status], GMatsReporter.__get_color(row[idx_status]))
            table_data.append(row)
        section.add_table(table_data, list(data_hits.columns), [('class', 'data')])
        section.add_paragraph('<b>Note:</b> Imperfect matches are marked with an asterisk (*) and are not used for determining the gMATS status.')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create overview table with status
        section = HtmlReportSection('gMATS')
        status = self._input_informs['gmats']['gMATS_status']
        section.add_table(
            [[HtmlTableCell(status, GMatsReporter.__get_color(status))]], ['gMATS status'], [('class', 'data')])

        # Create table with hits
        section.add_header('Detected hits', 3)
        self.__add_alleles_table(section)

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    @staticmethod
    def __get_color(status: str) -> str:
        """
        Returns the color for the corresponding gMATS status.
        :param status: gMATS status
        :return: Status cell
        """
        if status == 'covered':
            return 'green'
        if status == 'not_covered':
            return 'red'
        return 'yellow'
