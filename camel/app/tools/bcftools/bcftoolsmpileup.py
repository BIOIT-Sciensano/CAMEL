from pathlib import Path
from typing import Optional

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.toolpipeable import ToolPipeable


class BcftoolsMpileup(ToolPipeable):
    """
    Multi-way pileup producing genotype likelihoods.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools mpileup', '1.17')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("Reference genome input is required (FASTA)")
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError("Alignment input is required (BAM)")
        super()._check_input()

    def __get_output_key(self) -> str:
        """
        Returns the output key.
        :return: Output key
        """
        output_type = self._parameters['output_type'].value
        if output_type == 'b':
            return 'BCF_GZ'
        elif output_type == 'u':
            return 'BCF'
        elif output_type == 'z':
            return 'VCF_GZ'
        else:
            return 'VCF'

    def __get_output_path(self) -> Path:
        """
        Returns the path to the output file.
        :return: Output path
        """
        return self.folder / self._parameters['output_filename'].value

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command(self.__get_output_path())
        self._execute_command()
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(self.__get_output_path())]

    def _build_command(self, path_out: Optional[Path], pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command that is called.
        :param pipe_in: True if the tool is the input of a pipe
        :param pipe_out: True if the tool is the output of a pipe
        :param path_out: Output path
        :return: None
        """
        command_parts = [
            self._tool_command,
            str(self._tool_inputs['BAM'][0].path) if not pipe_in else '-',
            f"--fasta-ref {self._tool_inputs['FASTA'][0].path}"
        ]
        if 'BED_include' in self._tool_inputs:
            command_parts.append(f"--targets-file {self._tool_inputs['BED_include'][0].path}")
        elif 'BED_exclude' in self._tool_inputs:
            command_parts.append(f"--targets-file ^{self._tool_inputs['BED_exclude'][0].path}")
        command_parts.extend(self._build_options(excluded_parameters=['output_filename']))
        if not pipe_out:
            command_parts.extend(['--output', str(path_out)])
        self._command.command = ' '.join(command_parts)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _before_pipe(self, dir_: Path, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self._build_command(self.__get_output_path(), pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self._tool_outputs[self.__get_output_key()] = [ToolIOFile(self.__get_output_path())]
