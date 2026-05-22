import abc

from camel.app.core.command import Command
from camel.app.core.tool import Tool


class SeqkitBase(Tool, metaclass=abc.ABCMeta):
    """
    Baseclass for seqkit tools.
    """

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f"{self._tool_command.split(' ')[0]} version")
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()
