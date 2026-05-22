from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidParameterError, ToolExecutionError
from camel.app.tools.samtools.samtoolsbasepipeable import SamtoolsBasePipeable


class SamtoolsSort(SamtoolsBasePipeable):
    """
    Sorts alignment files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools sort', version=None)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(SamtoolsBasePipeable, self)._check_input()

    def _check_parameters(self) -> None:
        """
        Checks the tool parameters.
        :return: None
        """
        if self._parameters['output_format'].value.upper() not in ('SAM', 'BAM'):
            raise InvalidParameterError("Invalid output format (BAM/SAM supported)")
        super()._check_parameters()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self._check_stderr(self._command)

    def __build_command(self, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command
        :return: None
        """
        # Create excluded parameters
        excluded_params = ['output_filename'] if (pipe_out is True) else None

        # Construct command
        command_parts = [
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters=excluded_params))
        ]

        # Add input file
        if not pipe_in:
            command_parts.append(str(self._tool_inputs['BAM'][0].path))

        # Construct command
        self._command = Command(' '.join(command_parts))

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        output_path = self.folder / self._parameters['output_filename'].value
        if not output_path.is_file():
            raise ToolExecutionError(self.name, "Expected output not generated")
        output_key = self._parameters['output_format'].value.upper()
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully. Supersedes that of Tool class as samtools prints warnings to stderr.
        :return: None
        """
        if command.exit_code != 0:
            raise ToolExecutionError(self.name, f'Command execution failed (Exit code: {command.exit_code})')

    def _before_pipe(self, dir_, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self.__build_command(pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self.__set_output()
