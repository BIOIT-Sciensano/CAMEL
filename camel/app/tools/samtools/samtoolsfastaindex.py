from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import (
    InvalidParameterError,
    InvalidToolInputError,
    ToolExecutionError,
)
from camel.app.loggers import logger
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase


class SamtoolsFastaIndex(SamtoolsBase):
    """
    Indexes FASTA files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('samtools faidx', version=None)

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("No FASTA input file found")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError("Only one FASTA input file is supported.")
        super()._check_input()

    def _check_parameters(self) -> None:
        """
        Checks the parameters.
        :return: None
        """
        if 'regions' in self._parameters and 'output_filename' not in self._parameters:
            raise InvalidParameterError("Cannot extract regions without output filename")
        super()._check_parameters()

    def __symlink_input(self) -> Path:
        """
        Creates a symlink for the input.
        :return: Path to the symlink of the input
        """
        symlink_location = self.folder / self._tool_inputs['FASTA'][0].path.name
        expected_target = self._tool_inputs['FASTA'][0].path
        if symlink_location.is_symlink():
            if symlink_location.resolve() != expected_target.resolve():
                raise ToolExecutionError(
                    self.name,f"Symlink already exists but points to a different file: {symlink_location}")
            return symlink_location
        symlink_location.symlink_to(expected_target)
        return symlink_location

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if self.get_param_value('symlink_input') is True:
            fasta_in = self.__symlink_input()
        else:
            fasta_in = self._tool_inputs['FASTA'][0].path
        self.__build_command(fasta_in)
        self._execute_command()
        self._check_stderr(self._command)
        if 'regions' in self._parameters:
            self._tool_outputs['FASTA'] = [ToolIOFile(self.folder / self._parameters['output_filename'].value)]
        else:
            self._tool_outputs['FASTA'] = [ToolIOFile(fasta_in)]

    def __build_command(self, fasta_file: Path) -> None:
        """
        Builds the command for this tool.
        :param fasta_file: FASTA file
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, str(fasta_file)])
        if all(x in self._parameters for x in ('output_filename', 'regions')):
            logger.info("Extracting regions from FASTA file, file should already be indexed.")
            path_out = self._folder / self.get_param_value('output_filename')
            self._command.command += f" {self._parameters['regions'].value} > {path_out}"

    def _check_stderr(self, command: Command) -> None:
        """
        Checks the command stderr output.
        :param command: Command to check
        :return: None
        """
        if 'build FASTA index' in command.stderr:
            raise ToolExecutionError(self.name, "Cannot extract regions from an unindexed FASTA file.")
        super()._check_stderr(command)
