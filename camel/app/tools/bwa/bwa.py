import abc

from camelcore.app.command import Command

from camel.app.core.errors import ToolExecutionError
from camel.app.core.piping.toolpipeable import ToolPipeable


class BWA(ToolPipeable, metaclass=abc.ABCMeta):
    """
    Super class for reads mapping using BWA.
    """

    def _check_command_output(self, command: Command) -> None:
        """
        Parse stderr message of BWA cmd to check whether it runs successfully
        :return: None
        """
        if any(err in command.stderr.lower() for err in ('error', 'fail')):
            raise ToolExecutionError(self.name, f'Command failed: {command.command}\n stderr: {command.stderr}')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")
