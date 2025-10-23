from io import StringIO

import pandas as pd

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class Snpit(Tool):
    """
    Whole genome SNP-based identification of members of the Mycobacterium tuberculosis complex.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Snpit', '1.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['VCF'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_informs()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self.__generate_report_section())]

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '--input {}'.format(self._tool_inputs['VCF'][0].path)
        ])

    def __set_informs(self) -> None:
        """
        Sets the informs for this tool.
        :return: None
        """
        # noinspection PyTypeChecker
        data_all = pd.read_table(StringIO(self._command.stdout))
        data_output = data_all.to_dict('records')[0]
        self._informs['species'] = data_output['Species'] if not pd.isna(data_output['Species']) else 'NA'
        self._informs['lineage'] = data_output['Lineage'] if not pd.isna(data_output['Lineage']) else 'NA'
        self._informs['sublineage'] = data_output['Sublineage'] if not pd.isna(data_output['Sublineage']) else 'NA'
        self._informs['percent_matched'] = data_output['Percentage'] if not pd.isna(data_output['Percentage']) else 0.0

    def __generate_report_section(self) -> HtmlReportSection:
        """
        Creates a report section with the results of the tool.
        :return: Report section
        """
        section = HtmlReportSection('Snpit', 3)
        table_data = [
            ['Species:', '<i>{}</i>'.format(self._informs['species']) if self._informs['species'] != '' else 'NA'],
            ['Lineage:', self._informs['lineage'] if self._informs['lineage'] != '' else 'NA'],
            ['Sublineage:', self._informs['sublineage'] if self._informs['sublineage'] != '' else 'NA'],
            ['Percent match:', '{:.2f}%'.format(self._informs['percent_matched'])]
        ]
        section.add_table(table_data, table_attributes=[('class', 'information')])
        return section

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
