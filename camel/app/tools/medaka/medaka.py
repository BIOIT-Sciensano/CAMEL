import abc
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Medaka(Tool, metaclass=abc.ABCMeta):
    """
    Base class for Medaka tools.
    """

    def __init__(self, tool_name: str) -> None:
        """
        Initializes a Medaka tool.
        :return: None
        """
        super().__init__(tool_name, version=None)
        self._specific_parameters = []
        self._required_inputs = []
        self._input_string = ''
        self._output_string = ''
        self._option_string = ''
        self._output_type = ''

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f"{self._tool_command.split(' ')[0]} --version")
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _execute_tool(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        for key in self._required_inputs:
            if key not in self._tool_inputs:
                raise InvalidToolInputError(f"Input '{key}' is required")
        super()._check_input()

    def _set_input(self) -> None:
        """
        Sets the input specification.
        :return: None
        """
        pass

    def _set_output(self) -> None:
        """
        Sets the output specification.
        :return: None
        """
        self._output_string = str(self._parameters['output'])
        self._tool_outputs[self._output_type] = [ToolIOFile(Path(self._folder / self._parameters['output'].value))]

    def _set_specific_parameters(self) -> None:
        """
        Handles parameters that need specific handling for constructing the command line call.
        :return: None
        """
        pass

    def _build_command(self) -> None:
        """
        Build the command to run the tool.
        :return: None
        """
        self._option_string += " ".join(self._build_options(excluded_parameters=['output']))
        self._command.command = " ".join([
            self._tool_command, self._option_string, self._input_string, self._output_string
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_informs(self) -> None:
        """
        Sets the informs by analyzing the output.
        :return: None
        """
        pass
