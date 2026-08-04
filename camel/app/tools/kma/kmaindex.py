from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class KMAIndex(Tool):
    """
    Indexes a FASTA file for use with KMA ('kma index').
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('KMA index', None)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        super()._check_input()

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f"{self._tool_command.split()[0]} -v")
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split('-')[-1].strip()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self._tool_outputs['DB'] = [ToolIOValue(str(self.__output_prefix()))]

    def __output_prefix(self) -> Path:
        """
        Returns the path prefix of the generated KMA index.
        :return: Output prefix
        """
        return self.folder / self._tool_inputs['FASTA'][0].path.stem

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-i {self._tool_inputs['FASTA'][0].path}",
            f'-o {self.__output_prefix()}',
        ] + self._build_options())

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
