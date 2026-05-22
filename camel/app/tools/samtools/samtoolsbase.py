import abc

from camelcore.app.command import Command

from camel.app.core.errors import ToolExecutionError
from camel.app.core.tool import Tool


class SamtoolsBase(Tool, metaclass=abc.ABCMeta):
    """
    Super class for samtools.
    """

    def __init__(self, tool_name: str, version: str = None) -> None:
        """
        Initialize a samtools tool.
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
        command = Command(f'{self._tool_command.split()[0]} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split('\n')[0].split(' ')[-1].strip()

    def _check_stderr(self, command: Command) -> None:
        """
        Validate if the program ran correctly by checking the standard error.
        :param command: Command to check
        :return: None
        """
        if any(keyword in command.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError(self.name, command.stderr)
