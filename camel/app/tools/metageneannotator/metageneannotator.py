import os
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class MetaGeneAnnotator(Tool):
    """
    MetaGeneAnnotator is a gene-finding program for prokaryote and phage.
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('metageneannotator', '20080819')

    def _execute_tool(self):
        """
        Runs MetaGeneAnnotator
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
        super(MetaGeneAnnotator, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input keys given for MetaGeneAnnotator, '
                                                 'only FASTA allowed: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Only 1 input file allowed for MetaGeneAnnotator: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Only FASTA allowed as input key for MetaGeneAnnotator: {!r}'.format(self._tool_inputs))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        return os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)

    def __set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(Path(self.__get_basename() + '.out'))]

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self._tool_inputs['FASTA'][0].path
        options_string = ' '.join(self._build_options())
        outfile = self.__get_basename() + '.out'
        self._command.command = '{} {} {} > {}'.format(self._tool_command, input_string, options_string, outfile)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
