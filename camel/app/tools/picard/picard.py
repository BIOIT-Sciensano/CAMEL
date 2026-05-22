import abc
import re
from pathlib import Path

from camel.app.config import config
from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.piping.toolpipeable import ToolPipeable
from camel.app.loggers import logger


class Picard(ToolPipeable, metaclass=abc.ABCMeta):
    """
    Super class for Picard tools
    """

    def __init__(self, tool_name: str, version: str) -> None:
        """
        Initialize a picard tool
        :param tool_name: tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)

        self._function_name = None
        # parameters that should not be handled by self.build_options function
        self._specific_parameters = []
        self._required_inputs = ['BAM', 'SAM']
        self._main_input = ''
        self._input_string = ''
        self._output_string = ''
        self._output_type = 'BAM'
        # Elements for building command
        self._java_options = '-mx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false'
        self._java_options_temp_dir = f"TMP_DIR={config.dir_temp}"

    def update_java_options(self, java_options: str) -> None:
        """
        Returns the formatted java options of this tool.
        :return: Name
        """
        logger.info(f"Java options updated: '{java_options}'")
        self._java_options = f'{java_options}'

    def _execute_tool(self) -> None:
        """
        Function to run Picard function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        # BAM or SAM input required; mutually exclusive
        if ('SAM' in self._tool_inputs) and ('BAM' in self._tool_inputs):
            raise InvalidToolInputError()
        elif 'SAM' in self._tool_inputs:
            self._main_input = 'SAM'
            self._required_inputs.remove('BAM')
        elif 'BAM' in self._tool_inputs:
            self._main_input = 'BAM'
            self._required_inputs.remove('SAM')

        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidToolInputError(
                    f'Picard {self._name!r} required {input_file!r} input is missing in _tool_inputs!')
        super()._check_input()

    def _set_input(self) -> None:
        """
        Function to set main and extra inputs in self._input_string
        :return: None
        """
        if 'BAM' in self._tool_inputs:
            self._input_string += f"INPUT={self._tool_inputs['BAM'][0].path} "

        if 'SAM' in self._tool_inputs:
            self._input_string += f"INPUT={self._tool_inputs['SAM'][0].path} "

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f"R={self._tool_inputs['FASTA_REF'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles only one BAM file as output
        :return: None
        """
        self._tool_outputs[self._output_type] = [
            ToolIOFile(Path(self._folder) / self._parameters['output'].value)
        ]

    def _get_base_command(self) -> list[str]:
        """
        Returns the base command.
        :return: Base command
        """
        # Init base command
        if config.dependency_service == 'lmod':
            return ["java", self._java_options, "-jar $PICARD_JAR", self._tool_command, self._java_options_temp_dir]
        return ['picard', self._tool_command]

    def _build_command(self, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Build the command to run tool
        :param pipe_in: Tools receives input from pipe
        :param pipe_out: Tool outputs to pipe
        :return: None
        """
        command_parts = self._get_base_command()

        #Set input cmd line option
        if pipe_in:
            command_parts.append("I=/dev/stdin")
        else:
            command_parts.append(self._input_string)

        #Set output cmd line option
        if pipe_out:
            command_parts.append("O=/dev/stdout")
            self._specific_parameters.append('output')
        else:
            command_parts.append(self._output_string)

        #Add options
        command_parts.append(" ".join(self._build_options(excluded_parameters=self._specific_parameters, delimiter="=")))

        self._command.command = " ".join(command_parts)

    def _set_informs(self, stderr: str | None = None) -> None:
        """
        Analyse the result of picard run and update tool.informs, implement when necessary
        :return: None
        """
        pass

    def _check_command_output(self, command: Command) -> None:
        """
        Verify tool execution (return code) and stdout.
        :param command: Command to check
        :return: None
        """
        if command.exit_code != 0:
            raise ToolExecutionError(self.name, f'Picard {self._name,} fails to run, error msg: \n{command.stdout}')

        # log WARNINGs
        for line in command.stdout.splitlines():
            if re.match('WARNING', line):
                logger.warning(f' Picard - {line}')

    def _before_pipe(self, dir_: Path, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command(pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        self._set_informs(stderr)
