from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.tools.toolpipeable import ToolPipeable


class Star(ToolPipeable):
    """
    Super class for STAR indexing of reference genomes and alignment of spliced transcripts.
    """

    def __init__(self, tool_name: str, version: str) -> None:
        """
        Initializes a STAR tool.
        :param tool_name: Tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)
        self._required_inputs = []
        self._input_string = ""
        self._output_string = ""

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
        raise NotImplementedError("Method should be implemented by subclass.")

    def _build_command(self) -> None:
        """
        Builds the command to run STAR.
        :return: None
        """
        self._command.command = ' '.join([
          self._tool_command,
          self._input_string,
          *self._build_options(excluded_parameters=['filename_output', 'symlink_input']),
          self._output_string
        ])

    def _execute_tool(self) -> None:
        """
        Executes STAR.
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command()
        self._execute_command()

    def _set_output(self) -> None:
        """
        Sets the output specification and the output string.
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
