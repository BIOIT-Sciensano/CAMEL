from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Circos(Tool):
    """
    Circos is a software package for visualizing data and information. It visualizes data in a circular layout — this
    makes Circos ideal for exploring relationships between objects or positions.

    INPUT
        TXT: Text file with the configuration for Circos.
    OUTPUT
        PNG: Circos image rendered as PNG
        SVG: Circos image rendered as SVG
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Circos', '0.69-6')

    def _check_input(self) -> None:
        """
        Checks if the provided input is correct.
        :return: None
        """
        if 'TXT' not in self._tool_inputs:
            raise InvalidToolInputError("Circos configuration file ('TXT') is required.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, '-conf', str(self._tool_inputs['TXT'][0].path)])

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        png_path = self.folder / 'circos.png'
        if not png_path.is_file():
            raise ToolExecutionError(self.name, 'No PNG output file generated')
        self._tool_outputs['PNG'] = [ToolIOFile(png_path)]
        svg_path = self.folder / 'circos.svg'
        if not svg_path.is_file():
            raise ToolExecutionError(self.name, 'No SVG output file generated')
        self._tool_outputs['SVG'] = [ToolIOFile(svg_path)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
