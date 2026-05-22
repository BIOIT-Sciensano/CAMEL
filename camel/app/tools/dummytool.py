from camelcore.app.command import Command

from camel.app.core.tool import Tool
from camel.app.loggers import logger


class DummyTool(Tool):
    """
    Dummy tool for testing.
    """

    def __init__(self) -> None:
        """
        Initializes the dummy tool.
        :return: None
        """
        super().__init__('dummy', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes the dummy tool.
        """
        logger.info('Executing the dummy tool!')
        self._command = self.build_command()

    def build_command(self) -> Command:
        """
        Returns the command line call.
        :return: Command object
        """
        return Command(' '.join([
            self._tool_command,
            *self._build_options()
        ]))
