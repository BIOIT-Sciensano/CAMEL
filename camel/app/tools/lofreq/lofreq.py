import abc

from camel.app.core.command import Command
from camel.app.core.errors import ToolExecutionError
from camel.app.core.tool import Tool


class Lofreq(Tool, metaclass=abc.ABCMeta):
    """
    Super class for Lofreq.
    """

    def __init__(self, tool_name: str, version: str | None) -> None:
        """
        Initialize a Lofreq tool.
        :param tool_name: Tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{str(self._tool_command).split()[0]} version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split('\n')[0].split(':')[1].strip()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        super()._execute_tool()

    def _check_command_output(self, command: Command) -> None:
        """
        Validates if the program ran correctly by checking the standard error.
        :return: None
        """
        if self._command.stderr:
            if 'FATAL' in self._command.stderr:
                raise ToolExecutionError(self.name, f"{self.name} failed: '{self._command.stderr}'")
