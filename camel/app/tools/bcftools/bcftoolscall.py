from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidParameterError, InvalidToolInputError
from camel.app.core.piping.toolpipeable import ToolPipeable
from camel.app.loggers import logger
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsCall(BcftoolsBase, ToolPipeable):
    """
    SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools call', '1.17')

    def _check_parameters(self) -> None:
        """
        Checks the parameters.
        :return: None
        """
        if self._parameters['calling_method'].value not in ('consensus', 'multiallelic'):
            raise InvalidParameterError(f"Unrecognized snp calling method: {self._parameters['calling_method'].value}")
        if 'ploidy' not in self._parameters:
            logger.warning("Ploidy not specified will assume all sites are diploid.")

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('VCF', 'VCF_GZ', 'BCF')):
            raise InvalidToolInputError("No input file found (BCF / VCF / VCF_GZ supported)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command(self._get_output_path())
        self._execute_command()
        self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]

    def __get_input_file_path(self) -> Path:
        """
        Returns the path to the input file.
        :return: Input file path
        """
        return next(self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ', 'BCF') if k in self._tool_inputs)

    def _build_command(self, path_out: Path | None, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command.
        :param path_out: Output path
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        command_parts = [
            self._tool_command,
            '-' if pipe_in else str(self.__get_input_file_path()),
        ]

        # Calling method
        caller_flags = {
            'consensus': '--consensus-caller',
            'multiallelic': '--multiallelic-caller',
        }
        method = self.get_param_value('calling_method')
        if method not in caller_flags:
            raise ValueError(f"Unknown calling method: {method} (supported: {', '.join(caller_flags.keys())})")
        command_parts.append(caller_flags[method])

        if 'TXT_RG' in self._tool_inputs:
            command_parts.append(f'--regions-file {self._tool_inputs["TXT_RG"][0].path}')
        if 'TXT_SAMPLES' in self._tool_inputs:
            command_parts.append(f'--samples-file {self._tool_inputs["TXT_SAMPLES"][0].path}')

        command_parts += self._build_options(excluded_parameters=['calling_method', 'compress_output', 'output_filename'])
        if not pipe_out:
            command_parts.extend(['--output', str(path_out)])
        self._command.command = ' '.join(command_parts)

    def _before_pipe(self, dir_: Path, pipe_in: bool, pipe_out: bool) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self._build_command(self._get_output_path(), pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]
