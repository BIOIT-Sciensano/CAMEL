import os
import os.path
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class RefCov(Tool):
    """
    The RefCov software suite was written as a toolkit to provide multiple methods for analyzing coverage of sequence
    data across a reference.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('refcov', '0.3')

    def _execute_tool(self):
        """
        Runs RefCov
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - BAM key is required, BED key is optional
        - Only one input file allowed per key
        - No other input keys are allowed
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError(f'BAM input key is required for RefCov: {self._tool_inputs!r}')
        for key, value in self._tool_inputs.items():
            if key not in ['BAM', 'BED'] or len(value) > 1:
                raise InvalidToolInputError('Invalid input given for RefCov '
                                                     f'(BAM/BED and 1 file per key): {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError(f'Too many input keys given for RefCov (only BAM/BED allowed): {self._tool_inputs!r}')
        super()._check_input()

    def __build_basename(self, infiles):
        """
        Creates the prefix that will be used in the output
        :param infiles: List of input files
        :return: Prefix used in the output
        """
        return os.path.join(self._folder, '_'.join([os.path.splitext(f.path)[0] for f in infiles]))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs['BAM'][0].basename
        if 'BED' in self._tool_inputs:
            bedfile = self._tool_inputs['BED'][0].basename
            return self.__build_basename([infile, bedfile])
        else:
            return self.__build_basename([infile, infile])

    def __set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['TSV'] = [ToolIOFile(Path(basename + '_STATS.tsv'))]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        parts = ['--alignment-file-path {}'.format(self._tool_inputs['BAM'][0])]
        if 'BED' in self._tool_inputs:
            parts.append('--roi-file-path {}'.format(self._tool_inputs['BED'][0]))
        else:
            parts.append('--roi-file-path {}'.format(self._tool_inputs['BAM'][0]))
            parts.append('--roi-file-format bam')
        parts.append(f'--output-directory {self._folder}')
        return ' '.join(parts)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
