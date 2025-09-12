import pandas as pd

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class NanoPlot(Tool):
    """
    Plotting tool for long read sequencing data and alignments.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('NanoPlot', '1.41.6')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError('FASTQ input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        self._command.command = ' '.join([
            self._tool_command, '--fastq', str(self._tool_inputs['FASTQ'][0].path), *self._build_options()
        ])
        self._execute_command()
        self.__set_output()
        self.__set_informs()

    def __set_output(self) -> None:
        """
        Collects the tool output and stores it in tool outputs.
        """
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / 'NanoStats.txt')]
        self._tool_outputs['HTML'] = [ToolIOFile(self.folder / 'NanoPlot-report.html')]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks command output.
        :param command: Command to check
        :return: None
        """
        for line in command.stderr.splitlines():
            if 'skipping' in line:
                logger.warning(line)
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_informs(self) -> None:
        """
        Collects the informs.
        """
        data_summary = pd.read_table(self._tool_outputs['TSV'][0].path)
        for row in data_summary.to_dict('records'):
            self._informs[row['Metrics']] = row['dataset']
