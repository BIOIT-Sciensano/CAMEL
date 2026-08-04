from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Blat(Tool):
    """
    BLAT is a bioinformatics software tool which performs rapid mRNA/DNA and cross-species protein alignments. BLAT
    is more accurate and 500 times faster than popular existing tools for mRNA/DNA alignments and 50 times faster for
    protein alignments at sensitivity settings typically used when comparing vertebrate sequences.
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('blat', '36x2')

    def _execute_tool(self):
        """
        Runs BLAT
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - DB_DNA, DB_Prot or DB_DNAX required
        - FASTA_DNA, FASTA_Prot, FASTA_RNA, FASTA_DNAX, or FASTA_RNAX required
        - OOC optional
        - Only one input file allowed per input key
        :return: None
        """
        super()._check_input()
        if not 2 <= len(self._tool_inputs.keys()) <= 3:
            raise InvalidToolInputError(f'Invalid number (min 2, max 3) of input keys given for BLAT: {self._tool_inputs!r}')
        if not any(key in self._tool_inputs for key in ['DB_DNA', 'DB_Prot', 'DB_DNAX']):
            raise InvalidToolInputError(f'No database given for BLAT: {self._tool_inputs!r}')
        if not any(key in self._tool_inputs for key in ['FASTA_DNA', 'FASTA_Prot', 'FASTA_RNA', 'FASTA_DNAX', 'FASTA_RNAX']):
            raise InvalidToolInputError(f'No query file given for BLAT: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) == 3 and 'OOC' not in self._tool_inputs:
            raise InvalidToolInputError(f'Invalid input key given for BLAT: {self._tool_inputs!r}')
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidToolInputError(f'Too many input files per key (max = 1) given for BLAT: {self._tool_inputs!r}')

    def __get_input_key(self):
        """
        Returns the input key that is present in the inputs
        :return: Input key
        """
        return [x for x in self._tool_inputs if x in ['FASTA_DNA', 'FASTA_Prot', 'FASTA_RNA', 'FASTA_DNAX', 'FASTA_RNAX']][0]

    def __get_db_key(self):
        """
        Returns the database key key that is present in the inputs
        :return: Database key
        """
        return [x for x in self._tool_inputs if x in ['DB_DNA', 'DB_Prot', 'DB_DNAX']][0]

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs[self.__get_output_key()] = [
            ToolIOFile(self._folder / f'output{self.__get_output_extension()}')]

    def __get_output_key(self):
        """
        Returns the output key that has to be used based on the specified parameters in the database
        :return: Output key to use
        """
        output_key = 'TSV'
        if 'out' in self._parameters and self._parameters['out'].value not in ['psl', 'pslx', 'blast', 'blast8', 'blast9']:
            output_key = self._parameters['out'].value.upper()
        return output_key

    def __get_output_extension(self):
        """
        Returns the extension of the output file
        :return: Output file extension
        """
        return '.{}'.format(self._parameters['out'].value) if 'out' in self._parameters else '.psl'

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        return f'{self._tool_inputs[self.__get_db_key()][0]} {self._tool_inputs[self.__get_input_key()][0]}'

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options(delimiter=''))
        out_string = f'output{self.__get_output_extension()}'
        self._command.command = ' '.join([self._tool_command, input_string, options_string, out_string])

    def _check_command_output(self, command: Command):
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
