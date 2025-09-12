from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Filtlong(Tool):
    """
    Filtlong is a tool for filtering long reads by quality. It can take a set of long reads and produce a smaller,
    better subset. It uses both read length (longer is better) and read identity (higher is better) when choosing which
    reads pass the filter.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Filtlong', '0.2.0')

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
        # Construct command
        path_out = self.folder / FileSystemHelper.make_valid(self._parameters['output_name'].value)
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(excluded_parameters=['output_name']),
            str(self._tool_inputs['FASTQ'][0].path),
            f'> {path_out}'
        ])
        self._execute_command()

        # Collect output
        self._tool_outputs['FASTQ'] = [ToolIOFile(path_out)]
        self.__collect_stats(path_out)

    def __collect_stats(self, path_out: Path) -> None:
        """
        Collect statistics and stores them in the informs.
        :param path_out: Path to output file
        :return: None
        """
        self._informs['nb_reads_in'] = FastqUtils.count_reads(self._tool_inputs['FASTQ'][0].path)
        self._informs['nb_reads_out'] = FastqUtils.count_reads(path_out)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
