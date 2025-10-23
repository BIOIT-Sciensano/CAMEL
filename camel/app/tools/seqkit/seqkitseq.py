from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils, fastautils, fastqutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.seqkit.seqkitbase import SeqkitBase


class SeqkitSeq(SeqkitBase):
    """
    Seqkit seq performs common transformations of FASTA / FASTQ files.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self) -> None:
        """
        Initialize seqkit seq.
        :return: None
        """
        super().__init__(name='Seqkit seq', version=None)

    def __get_input_key(self) -> str:
        """
        Returns the key for the input file.
        :return: Input key
        """
        return 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in SeqkitSeq.INPUT_KEYS):
            raise InvalidToolInputError('{} input is required.'.format(' or '.join(SeqkitSeq.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self._folder / self._parameters['output_filename'].value
        self.__build_command()
        self._execute_command()
        self.__set_output(output_path)
        self.__collect_stats(output_path)
        self._informs['min_length'] = self._parameters['min_length'].value if 'min_length' in self._parameters else None
        self._informs['min_qual'] = self._parameters['min_qual'].value if 'min_qual' in self._parameters else None

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        input_key = self.__get_input_key()
        self._command.command = ' '.join([
            self._tool_command,
            *[str(f.path) for f in self._tool_inputs[input_key]],
            *self._build_options()
        ])

    def __set_output(self, output_path: Path) -> None:
        """
        Sets the tool output.
        :param output_path: Output path
        :return: None
        """
        output_key = 'FASTA' if output_path.name.lower().endswith('.fasta') else 'FASTQ'
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __collect_stats(self, path_out: Path) -> None:
        """
        Collect statistics and stores them in the informs.
        :param path_out: Path to output file
        :return: None
        """
        if self.__get_input_key() == 'FASTQ':
            self._informs['nb_seqs_in'] = fastqutils.count_reads(self._tool_inputs['FASTQ'][0].path)
            self._informs['nb_seqs_out'] = fastqutils.count_reads(path_out)
        else:
            self._informs['nb_seqs_in'] = fastautils.count_reads(self._tool_inputs['FASTA'][0].path)
            self._informs['nb_seqs_out'] = fastautils.count_reads(path_out)
