from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class FragGeneScan(Tool):
    """
    FragGeneScan is an application for finding (fragmented) genes in short reads. It can also be applied to predict
    prokaryotic genes in incomplete assemblies or complete genomes.
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('fraggenescan', '1.30')

    def _execute_tool(self):
        """
        Runs FragGeneScan
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only one input file allowed per key
        - No other input keys are allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs or len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError(f'Invalid input keys given for FragGeneScan, only FASTA allowed: {self._tool_inputs!r}')
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError(f'Invalid number of files given for FragGeneScan, only 1 allowed: {self._tool_inputs!r}')

    def __get_basename(self) -> Path:
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        return Path(self._folder, self._tool_inputs['FASTA'][0].basename)

    def __set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['TSV'] = [ToolIOFile(Path(str(basename) + '.processed.out'))]
        self._tool_outputs['FASTA'] = [ToolIOFile(Path(str(basename) + '.processed.ffn'))]
        self._tool_outputs['FASTA_Prot'] = [ToolIOFile(Path(str(basename) + '.processed.faa'))]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        parts = ['-genome={}'.format(self._tool_inputs['FASTA'][0]),
                 '-out={}'.format(str(self.__get_basename()) + '.processed')]
        return ' '.join(parts)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options(delimiter='='))
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
