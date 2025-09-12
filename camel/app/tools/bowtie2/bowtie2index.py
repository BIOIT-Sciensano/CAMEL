from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.components.files.fileutils import FileUtils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class Bowtie2Index(Tool):
    """
    Index genome using 'bowtie2-build' cmd of Bowtie2
    """

    MULTI_FASTA_GENOME_FILE = 'concatenated.fasta'

    def __init__(self) -> None:
        """
        Initialize bowtie2 index
        :return: None
        """
        super().__init__('bowtie2 index', '2.5.1')

    def _execute_tool(self) -> None:
        """
        Function to run BWA index
        :return: None
        """
        path_fasta = self.__collect_fasta_input()
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            str(path_fasta),
            str(path_fasta)])
        self._execute_command()
        self._tool_outputs['INDEX_GENOME_PREFIX'] = [ToolIOValue(str(path_fasta))]

    def __get_multi_fasta_genome_filename(self) -> Path:
        """
        Get the filename used for multi fasta file representing complete genome
        :return: name of the multi fasta file with complete path
        """
        return self._folder / Bowtie2Index.MULTI_FASTA_GENOME_FILE

    def _check_input(self) -> None:
        """
        Check FASTA_REF input and concatenate them if multiple fasta input files
        :return: None
        """
        if 'FASTA_REF' not in self._tool_inputs or len(self._tool_inputs['FASTA_REF']) == 0:
            raise InvalidToolInputError("FASTA_REF input is required")
        super()._check_input()

    def __collect_fasta_input(self) -> Path:
        """
        Collects the FASTA input
        :return: Path to FASTA file that will be indexed.
        """
        nb_of_inputs = len(self._tool_inputs['FASTA_REF'])
        if nb_of_inputs > 1:
            logger.info(f'Creating concatenated FASTA file ({nb_of_inputs} files)')
            path_multi_fasta = self.folder / Bowtie2Index.MULTI_FASTA_GENOME_FILE
            FileUtils.concatenate_files(path_multi_fasta, [f.path for f in self._tool_inputs['FASTA_REF']])
            return path_multi_fasta
        else:
            path_fasta = self._folder / self._tool_inputs['FASTA_REF'][0].path.name
            if not path_fasta.exists():
                logger.info(f'Creating symlink for input FASTA file: {path_fasta}')
                path_fasta.symlink_to(self._tool_inputs['FASTA_REF'][0].path)
            return path_fasta

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully.
        :param command: Command to execute
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
