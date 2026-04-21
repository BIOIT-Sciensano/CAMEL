import json
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils, fileutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class SeqSero2(Tool):
    """
    Salmonella serotype prediction from genome sequencing data.
    """

    def __init__(self) -> None:
        """
        Initialize tool.
                :return: None
        """
        super().__init__('SeqSero2', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _execute_tool(self) -> None:
        """
        Execute the tool.
        :return: None
        """
        self.build_command(self._parameters['mode'].value)
        self._execute_command()
        self.__set_output()
        input_folder = self._tool_inputs['DIR'][0].path
        self._informs['_tag'] = self._parameters['mode'].value
        self.__add_informs(input_folder)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super()._check_input()
        if not self._parameters.get('mode') \
                or self._parameters['mode'].value not in ('kmer', 'allele', 'kmerread'):
            raise InvalidToolInputError("A Seqsero2 processing mode must be passed to the tool, "
                                                 "choose from kmer, allele, or kmerread")
        if self._parameters['mode'].value == 'kmer':
            if 'FASTA' not in self._tool_inputs:
                raise InvalidToolInputError("FASTA input is required in kmer mode")
        else:
            if sum(x in self._tool_inputs for x in ('FASTQ', 'FASTQ_PE', 'FASTQ_ONT')) != 1:  # not exactly one
                raise InvalidToolInputError(f"Exactly one FASTQ input is required in "
                                                     f"{self._parameters['mode'].value} mode")
            if self._parameters['mode'].value == 'allele' and 'FASTQ_ONT' in self._tool_inputs:
                raise InvalidToolInputError("allele mode is not available for nanopore FQ input.")

    def __set_output(self) -> None:
        """
        Collects the tool output files.
        :return: None
        """
        self._tool_outputs['TXT'] = [ToolIOFile(self.folder / 'SeqSero_result.txt')]

    def build_command(self, mode: str) -> None:
        """
        Concatenates required parameters and options to build the command.
        :param mode: seqsero2 execution mode; 'allele' or 'kmer' or 'kmerread'
        :return: None
        """
        command_parts = [self._tool_command, '-d', str(self.folder), " ".join(self._build_options(excluded_parameters=['mode']))]
        if mode == 'kmer':
            command_parts.extend(['-t 4 -m k', '-i', str(self._tool_inputs['FASTA'][0])])
        else:
            if mode == 'allele':
                command_parts.append('-m a')
            else:  # if mode == 'kmerread':
                command_parts.append('-m k')

            if 'FASTQ_ONT' in self._tool_inputs:  # in case of ONT input data
                # create intermediary input dir because Seqsero2 needs a different input than output dir
                (self.folder / 'in').mkdir()
                if fileutils.is_gzipped(self._tool_inputs['FASTQ_ONT'][0].path):
                    # Gunzip file because in -t 5 the input needs to be gunzipped.
                    fastq_gunzipped = self.folder / 'in' / self._tool_inputs['FASTQ_ONT'][0].path.stem
                    fileutils.gzip_extract(self._tool_inputs['FASTQ_ONT'][0].path, fastq_gunzipped)
                    command_parts.extend(['-t 5', '-i', str(fastq_gunzipped)])
                else:
                    command_parts.extend(['-t 5', '-i', str(self._tool_inputs['FASTQ_ONT'][0].path)])
            else:  # if 'FASTQ_PE' in self._tool_inputs:
                command_parts.extend([
                    '-t 2', '-i', str(self._tool_inputs['FASTQ_PE'][0].path), str(self._tool_inputs['FASTQ_PE'][1].path)])
        self._command.command = ' '.join(command_parts)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __add_informs(self, input_folder: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param input_folder: Input database directory
        :return: None
        """
        path_metadata = input_folder / 'db_update_info.json'
        if not path_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {path_metadata}')
        with path_metadata.open() as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        self._informs['db_path'] = str(input_folder)
