import abc

from camelcore.app.command import Command

from camel.app.core.tool import Tool


class Minimap2Base(Tool, metaclass=abc.ABCMeta):
    """
    Baseclass for minimap2 tools.
    """

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.strip().split('-')[0]
