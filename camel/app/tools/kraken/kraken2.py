from datetime import datetime
from pathlib import Path

import re

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class Kraken2(Tool):
    """
    Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic
    studies. Previous attempts by other bioinformatics software to accomplish this task have often used sequence
    alignment or machine learning techniques that were quite slow, leading to the development of less sensitive but
    much faster abundance estimation programs. Kraken aims to achieve high sensitivity and high speed by utilizing
    exact alignments of k-mers and a novel classification algorithm.
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('kraken2', '2.1.1')
        self._input_key = None

    def _execute_tool(self) -> None:
        """
        Runs Kraken 2.
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__set_informs(self._command.stderr)

    def _check_input(self) -> None:
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
                raise InvalidToolInputError(
                    f'There is more than 1 FASTA/Q file or more/less than two FASTQ_PE files given for Kraken: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError(
                f'Too many input keys given for Kraken ((FASTA or FASTQ or FASTQ_PE) and DB): {self._tool_inputs!r}')

    def __get_basename(self) -> str:
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = Path(self._tool_inputs[self._input_key][0].path)
        return str(Path(self._folder) / infile.stem)

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        basename = Path(self.__get_basename())
        self._tool_outputs['TSV'] = [ToolIOFile(basename.parent / f'{basename.name}.output.tsv')]
        self._tool_outputs['TSV_report'] = [ToolIOFile(basename.parent / f'{basename.name}.report.tsv')]

    def __build_input_string(self) -> str:
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        command_parts = ['--db', str(self._tool_inputs['DB'][0].path)]
        if self._input_key == 'FASTA':
            command_parts.append(str(self._tool_inputs['FASTA'][0].path))
        elif self._input_key == 'FASTQ':
            command_parts.append(str(self._tool_inputs['FASTQ'][0].path))
        else:
            command_parts.extend(
                [str(self._tool_inputs['FASTQ_PE'][0].path), str(self._tool_inputs['FASTQ_PE'][1].path), '--paired'])
        command_parts.extend([
            '--output {}'.format(self.__get_basename() + '.output.tsv'),
            '--report {}'.format(self.__get_basename() + '.report.tsv')
        ])
        return ' '.join(command_parts)

    def __set_input_key(self) -> None:
        """
        Sets the instance variable self._input_key
        :return: None
        """
        for key in self._tool_inputs.keys():
            if key != 'DB':
                self._input_key = key

    def __set_informs(self, stderr: str) -> None:
        """
        Sets the informs based on the input database.
        :return: None
        """
        for line in stderr.splitlines():
            m = re.match(r'.*Using database (.*) created with RefSeq sequences from (\d+)', line.strip())
            if m:
                date_update = datetime.strptime(m.group(2), '%Y%m%d').strftime('%d-%m-%Y')
                self._informs['database'] = {'name': m.group(1), 'last_update': date_update}
                return
        self._informs['database'] = {'name': 'n/a', 'last_update': 'n/a'}

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        options_string = ' '.join(self._build_options())
        self._command.command = f'{self._tool_command} {self.__build_input_string()} {options_string}'

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        toolutils.check_tool_execution(self, command, exit_code=0)
