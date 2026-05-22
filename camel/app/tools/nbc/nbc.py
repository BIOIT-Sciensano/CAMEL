import gzip
import os
import tarfile
from typing import IO

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils


class Nbc(Tool):
    """
    The naive Bayes classifier (NBC) classifies metagenomic reads to their best taxonomic match. Results indicate that
    NBC can assign next- generation sequencing reads to their taxonomic classification and can find significant
    populations of genera that other classifiers may miss.
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('NBC', '1.1')

    def _execute_tool(self) -> None:
        """
        Runs NBC
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__tabulate_output()
        self.__process_output(self.__get_read_lengths())
        self.__archive_raw_results()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - Only FASTA and DB are allowed and required
        - Only one input file allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs or 'DB' not in self._tool_inputs:
            raise InvalidToolInputError(
                f'Invalid input keys given for NBC, FASTA and DB required: {self._tool_inputs!r}'
            )
        if len(self._tool_inputs.keys()) != 2:
            raise InvalidToolInputError(
                'Invalid number of input keys given for NBC, '
                f'only FASTA and DB allowed: {self._tool_inputs!r}'
            )
        if len(self._tool_inputs['FASTA']) > 1 or len(self._tool_inputs['DB']) > 1:
            raise InvalidToolInputError(
                f'Invalid number (max = 1) of files per key given for NBC: {self._tool_inputs!r}'
            )

    def __set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['TGZ'] = [ToolIOFile(self._folder / 'raw_output.tar.gz')]
        self._tool_outputs['TSV'] = [ToolIOFile(self._folder / 'processed_output.tsv')]

    def __build_input_string(self) -> str:
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        parts = [
            f'-a {self._tool_inputs["FASTA"][0]}',
            f'-j {self._tool_inputs["DB"][0]}',
        ]
        return ' '.join(parts)

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join(
            [self._tool_command, input_string, options_string]
        )

    def __tabulate_output(self) -> None:
        """
        Creates the tabulated output from the raw output files
        :return: None
        """
        self.__create_dir_list()
        self._command.command = 'tabulate -g genomes.txt *.txt.gz'
        self._execute_command()

    def __create_dir_list(self) -> None:
        """
        Creates a file containing all directories in the database
        :return: None
        """
        dir_list = os.listdir(self._tool_inputs['DB'][0].path)
        with open(
            os.path.join(self._folder, 'genomes.txt'), 'w', encoding='utf-8'
        ) as outfile:
            for entry in dir_list:
                outfile.write(entry + '\n')

    def __process_output(self, read_lengths: dict[str, int]) -> None:
        """
        Converts the raw output to a table with more information (i.e. the read name, the raw score, the significance
        cutoff, whether the score is significant (True/False) and the species name that was identified
        :param read_lengths: Dictionary mapping read names to their lengths
        :return: None
        """
        # Example output file
        # names Read1        Read2        Read3        Read4        Read5
        # name1 -4394.422894 -4481.248143 -4466.777268 -4449.821487 -4481.248143
        # name2 -4479.885575 -4419.27694  -4449.581257 -4495.037734 -4495.037734
        # name3 -4495.461521 -4495.461521 -4495.461521 -4466.004776 -4466.004776
        with open(
            os.path.join(self._folder, 'processed_output.tsv'), 'w', encoding='utf-8'
        ) as outf:
            outf.write(
                '\t'.join(['Read name', 'Score', 'Cutoff', 'Significant', 'Species\n'])
            )
            for entry in os.listdir(self._folder):
                if entry.endswith('.csv.gz'):
                    with gzip.open(os.path.join(self._folder, entry), 'rb') as infile:
                        lines = []
                        for line in infile:
                            lines.append(self.__to_floats(line))
                        zipped_lines = list(zip(*lines))
                        names = zipped_lines[0]
                        for line in zipped_lines[1:]:
                            # Find the entry with the highest value as this is the best hit
                            index = line.index(max(line[1:]))
                            self.__write_to_output(
                                outf, line, names, index, read_lengths
                            )

    @staticmethod
    def __to_floats(line: bytes) -> list[float | str]:
        """
        Creates a list of the columns in the given line with all numeric values converted to float (i.e. all columns
        except the name column and the first row).
        :param line: Line to split
        :return: List of column values
        """
        items = line.replace(b'-inf', b'-999999999999999.0').strip().split(b',')
        new_items: list[float | str] = []
        for item in items:
            try:
                new_items.append(float(item))
            except ValueError:
                new_items.append(item.decode())
        return new_items

    def __write_to_output(
        self,
        outf: IO[str],
        line: tuple[float | str, ...],
        names: tuple[float | str, ...],
        index: int,
        read_lengths: dict[str, int],
    ) -> None:
        """
        Writes the final output based on the given parameters
        :param outf: Output file
        :param line: Line from the raw NBC tabulated output file
        :param names: List of read names
        :param index: Index of the maximum value
        :param read_lengths: Lengths of the reads to calculate significance
        :return: None
        """
        result = [
            line[0],
            str(line[index]),
            str(self.__get_cutoff(read_lengths[line[0]])),
            self.__is_significant(line[index], read_lengths[line[0]]),
            names[index],
        ]
        outf.write('\t'.join(result))
        outf.write('\n')

    def __is_significant(self, value: float, length: int) -> str:
        """
        Returns whether a value can be considered significant (i.e. the species is found)
        :param value: Value to check
        :param length: Length of the read
        :return: 'True' or 'False'
        """
        return 'True' if value > self.__get_cutoff(length) else 'False'

    @staticmethod
    def __get_cutoff(length: int) -> float:
        """
        Returns the cutoff that can be used to see whether an assignment is significant
        :param length: Length of the read
        :return: Cutoff
        """
        return -23.7 * length + 490

    def __get_read_lengths(self) -> dict[str, int]:
        """
        Stores all read lengths in a dictionary with the read name as key and the length as value
        :return: Dictionary of read lengths
        """
        with open(self._tool_inputs['FASTA'][0].path, encoding='utf-8') as infile:
            read_lengths: dict[str, int] = {}
            current_id = ''
            for line in infile:
                if line.startswith('>'):
                    current_id = line[1:].strip()
                    read_lengths[current_id] = 0
                else:
                    read_lengths[current_id] += len(line.strip())
        return read_lengths

    def __archive_raw_results(self) -> None:
        """
        Archives the raw output files in a tar.gz file
        :return: None
        """
        with tarfile.open(
            os.path.join(self._folder, 'raw_output.tar.gz'), 'w:gz'
        ) as tar:
            for entry in os.listdir(self._folder):
                if entry.endswith('.csv.gz'):
                    gzfile = os.path.join(self._folder, entry)
                    with (
                        gzip.open(gzfile, 'rb') as gz,
                        open(gzfile[:-3], 'wb') as tempf,
                    ):
                        for line in gz:
                            tempf.write(line)
                    tar.add(gzfile[:-3], arcname=os.path.basename(gzfile[:-3]))
                    os.remove(gzfile[:-3])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
