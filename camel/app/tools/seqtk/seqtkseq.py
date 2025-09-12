from pathlib import Path

from camel.app.components.files.fastautils import FastaUtils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class SeqtkSeq(Tool):
    """
    Seqtk seq perform common transformations of FASTA / FASTQ files.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self) -> None:
        """
        Initialize seqtk seq.
        :return: None
        """
        super().__init__('Seqtk seq', '1.4')

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
        if not any(x in self._tool_inputs for x in SeqtkSeq.INPUT_KEYS):
            raise InvalidToolInputError('{} input is required.'.format(' or '.join(SeqtkSeq.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self._folder / self._parameters['output_filename'].value
        self.__build_command(output_path)
        self._execute_command()
        self.__set_output(output_path)
        self.__collect_stats(output_path)

    def __build_command(self, output_path: Path) -> None:
        """
        Builds the command line call.
        :param output_path: Output file path
        :return: None
        """
        input_key = self.__get_input_key()
        self._command.command = ' '.join([
            self._tool_command,
            *[str(f.path) for f in self._tool_inputs[input_key]],
            *self._build_options(['output_filename']),
            f'> {output_path}'
        ])

    def __set_output(self, output_path: Path) -> None:
        """
        Sets the tool output.
        :param output_path: Output path
        :return: None
        """
        output_key = 'FASTA' if output_path.name.lower().endswith('.fasta') else 'FASTQ'
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def __collect_stats(self, path_out: Path) -> None:
        """
        Collect statistics and stores them in the informs.
        :param path_out: Path to output file
        :return: None
        """
        if self.__get_input_key() == 'FASTQ':
            self._informs['nb_seqs_in'] = FastqUtils.count_reads(self._tool_inputs['FASTQ'][0].path)
            self._informs['nb_seqs_out'] = FastqUtils.count_reads(path_out)
        else:
            self._informs['nb_seqs_in'] = FastaUtils.count_reads(self._tool_inputs['FASTA'][0].path)
            self._informs['nb_seqs_out'] = FastaUtils.count_reads(path_out)
