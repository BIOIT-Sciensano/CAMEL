import abc

from camel.app.command.command import Command
from camel.app.error import InvalidToolInputError, ToolExecutionError
from camel.app.tools.tool import Tool


class Bedtools(Tool, metaclass=abc.ABCMeta):
    """
    The base class for bedtools.
    """

    def __init__(self, tool_name: str, version: str) -> None:
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)
        self._required_inputs = []

    def _check_command_output(self, command: Command) -> None:
        """
        Validate if the program ran correctly by checking the standard error.
        :param command: Command
        :return: None
        """
        if any(keyword in command.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError(self.name, f'{self.name} failed: {command.stderr}')

    def _check_required_inputs(self) -> None:
        """
        Check required input.
        :return: None
        """
        if len(self._required_inputs) == 0:
            return
        try:
            next(key for key in self._required_inputs if key in self._tool_inputs)
        except StopIteration:
            raise InvalidToolInputError(f"Required inputs missing ({', '.join(self._required_inputs)})")
