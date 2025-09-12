import abc
from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidParameterError
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class BcftoolsBase(Tool, metaclass=abc.ABCMeta):
    """
    Baseclass for bcftools.
    """

    def _get_output_key(self) -> str:
        """
        Returns the key for the output file.
        :return: Output file key
        """
        if 'output_type' not in self._parameters:
            return 'VCF'
        try:
            return {
                'b': 'BCF_GZ',
                'u': 'BCF',
                'z': 'VCF_GZ',
                'v': 'VCF'
            }[self._parameters['output_type'].value]
        except KeyError:
            raise InvalidParameterError(f"Invalid output format: {self._parameters['VCF'].value}")

    def _get_output_path(self) -> Path:
        """
        Returns the output filename.
        :return: Output filename
        """
        if 'output_filename' in self._parameters:
            return self.folder / str(self._parameters['output_filename'].value)
        logger.info("Output filename not set, reverting to default.")
        return self.folder / f"bcftools_out.{self._get_output_key().lower().replace('_', '.')}"

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully
        :param command: Command to check.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
