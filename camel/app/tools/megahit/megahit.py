import os

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.tools.tool import Tool
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile


class Megahit(Tool):
    """
    MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('megahit', '1.1.1-2')
        self.__input_key = None

    def _execute_tool(self):
        """
        Runs Megahit
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA or FASTQ allowed either SE, PE, or INT
        - One input file allowed for SE or INT, two files for PE
        :return: None
        """
        super()._check_input()
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Invalid number of input keys given for Megahit '
                                                 '(only 1 allowed): {!r}'.format(self._tool_inputs))
        if list(self._tool_inputs.keys())[0] not in ['FASTQ_PE', 'FASTA_PE', 'FASTQ_SE', 'FASTA_SE', 'FASTQ_INT', 'FASTA_INT']:
            raise InvalidToolInputError('Not enough valid input files given for Megahit '
                                                 '(only FASTQ/A - SE, PE, or INT allowed): {!r}'.format(self._tool_inputs))
        key, value = list(self._tool_inputs.items())[0]
        if (key.endswith('PE') and len(value) != 2) or (not key.endswith('PE') and len(value) != 1):
            raise ValueError('Invalid number (2 for PE, 1 for SE or INT) of files per key given '
                             'for Megahit: {!r}'.format(self._tool_inputs))

    def __set_input_key(self):
        """
        Sets the instance variable self.__input_key
        :return: None
        """
        self.__input_key = list(self._tool_inputs.keys())[0]

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / 'output' / 'final.contigs.fa')]
        self._tool_outputs['LOG'] = [ToolIOFile(self._folder/ 'output' / 'log')]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        inputs = {
            'PE': lambda: '-1 {} -2 {}'.format(self._tool_inputs[self.__input_key][0],
                                               self._tool_inputs[self.__input_key][1]),
            'SE': lambda: '--read {}'.format(self._tool_inputs[self.__input_key][0]),
            'INT': lambda: '--12 {}'.format(self._tool_inputs[self.__input_key][0])
        }
        return inputs[self.__input_key.split('_')[1]]()

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        out_string = '--out-dir {}'.format(os.path.join(self._folder, 'output'))
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, out_string, options_string])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
