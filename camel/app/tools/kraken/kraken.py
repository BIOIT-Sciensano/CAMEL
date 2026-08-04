import os
import re
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Kraken(Tool):
    """
    Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic
    studies. Previous attempts by other bioinformatics software to accomplish this task have often used sequence
    alignment or machine learning techniques that were quite slow, leading to the development of less sensitive but
    much faster abundance estimation programs. Kraken aims to achieve high sensitivity and high speed by utilizing
    exact alignments of k-mers and a novel classification algorithm.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('kraken', '0.10.5')
        self._input_key = None

    def _execute_tool(self):
        """
        Runs Prinseq
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - (FASTA or FASTQ or FASTQ_PE) and DB keys are required
        - Only one input file allowed per key (2 for FASTQ_PE)
        - No other input keys are allowed
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('FASTA', 'FASTQ', 'FASTQ_PE')) or 'DB' not in self._tool_inputs:
            raise InvalidToolInputError(f'FASTA/Q input or DB input missing for Kraken: {self._tool_inputs!r}')
        for key, value in self._tool_inputs.items():
            if (key != 'FASTQ_PE' and len(value) > 1) or (key == 'FASTQ_PE' and len(value) != 2):
                raise InvalidToolInputError('There is more than 1 FASTA/Q file or more/less than two FASTQ_PE '
                                                     f'files given for Kraken: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError('Too many input keys given for Kraken ((FASTA or FASTQ or FASTQ_PE) '
                                                 f'and DB): {self._tool_inputs!r}')

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs[self._input_key][0].basename
        return os.path.join(self._folder, os.path.splitext(infile)[0])

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['TSV'] = [ToolIOFile(Path(basename + '.output.tsv'))]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        command_parts = []
        if self._input_key == 'FASTA':
            command_parts.append('{}'.format(self._tool_inputs['FASTA'][0]))
        elif self._input_key == 'FASTQ':
            command_parts.append('{} {}'.format('--fastq-input', self._tool_inputs['FASTQ'][0]))
        else:
            command_parts.append('{} {} {}'.format(self._tool_inputs['FASTQ_PE'][0], self._tool_inputs['FASTQ_PE'][1],
                                                   '--paired --fastq-input'))
        command_parts += ['--db {}'.format(self._tool_inputs['DB'][0]),
                          '--output {}'.format(self.__get_basename() + '.output.tsv')]
        return ' '.join(command_parts)

    def __set_input_key(self):
        """
        Sets the instance variable self._input_key
        :return: None
        """
        for key in self._tool_inputs.keys():
            if key != 'DB':
                self._input_key = key

    def __check_preload(self, options):
        """
        Checks whether the provided database is located at '/dev/shm' or not and adds the preload option accordingly
        :param options: Option string to check
        :return: Option string
        """
        if self._tool_inputs['DB'][0].path.name.startswith('/dev/shm'):
            return re.sub(r'--preload\s', '', options)
        elif 'preload' not in options:
            return '--preload ' + options
        raise ValueError('Cannot determine pre-load')

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        options_string = self.__check_preload(' '.join(self._build_options()))
        self._command.command = f'{self._tool_command} {self.__build_input_string()} {options_string}'

    def _check_command_output(self, command: Command):
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        if command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Command execution failed (Exit code: {command.exit_code})")
