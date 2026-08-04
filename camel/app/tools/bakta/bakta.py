from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Bakta(Tool):
    """
    Bakta is a tool for the rapid & standardized annotation of bacterial genomes
    and plasmids from both isolates and MAGs.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Bakta', '1.9.4')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError("FASTA input requires exactly 1 file.")
        super()._check_input()

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['FASTA'][0].path),
            *self._build_options()]
        )

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self) -> None:
        """
        Collects the output files of interest.
        :return: None
        """
        fasta_file = self._tool_inputs['FASTA'][0]
        filename = fasta_file.path.stem
        self._tool_outputs['FAA'] = [ToolIOFile(self.folder / self._parameters['output_dir'].value / f'{filename}.faa')]
        self._tool_outputs['GFF3'] = [ToolIOFile(self.folder / self._parameters['output_dir'].value / f'{filename}.gff3')]
        self._tool_outputs['GBFF'] = [ToolIOFile(self.folder / self._parameters['output_dir'].value / f'{filename}.gbff')]
