from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class PharokkaMultiplotter(Tool):
    """
    Pharokka Multiplotter is a companion functionality of Pharokka to visualize genome annotations.
    """

    OUTPUT_DICT = {
        'PNG': '.png',
        'SVG': '.svg'
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Pharokka', '1.7.3')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()

        # Collect the output
        path_out = Path(self.folder, self._parameters['outdir'].value)
        for key, extension in PharokkaMultiplotter.OUTPUT_DICT.items():
            for path_plot in sorted(path_out.glob(f'*{extension}')):
                if key not in self._tool_outputs:
                    self._tool_outputs[key] = []
                self.tool_outputs[key].append(ToolIOFile(path_plot))

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'GBK' not in self._tool_inputs:
            raise InvalidToolInputError("No genbank input found")
        if len(self._tool_inputs['GBK']) != 1:
            raise InvalidToolInputError("Only one genbank file can be annotated at a time.")
        super()._check_input()

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __build_command(self) -> None:
        """
        Build the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            '-g',
            str(self._tool_inputs['GBK'][0].path)
        ])
