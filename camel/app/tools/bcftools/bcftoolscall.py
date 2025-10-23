from pathlib import Path

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import InvalidParameterError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsCall(BcftoolsBase):
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
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_input_file_path(self) -> Path:
        """
        Returns the path to the input file.
        :return: Input file path
        """
        return next(self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ', 'BCF') if k in self._tool_inputs)

    def __build_command(self) -> None:
        """
        Builds the command.
        :return: None
        """
        command_parts = [
            self._tool_command,
            str(self.__get_input_file_path()),
        ]
        if self._parameters['calling_method'].value == 'consensus':
            command_parts.append('--consensus-caller')
        else:
            command_parts.append('--multiallelic-caller')

        if 'TXT_RG' in self._tool_inputs:
            command_parts.append(f'--regions-file {self._tool_inputs["TXT_RG"][0].path}')
        if 'TXT_SAMPLES' in self._tool_inputs:
            command_parts.append(f'--samples-file {self._tool_inputs["TXT_SAMPLES"][0].path}')

        command_parts += self._build_options(['calling_method', 'compress_output'])
        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]
