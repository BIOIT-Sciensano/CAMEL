import os
import os.path
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.tool import Tool
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile


class RefCovSort(Tool):
    """
    The RefCov software suite was written as a toolkit to provide multiple methods for analyzing coverage of sequence
    data across a reference. This Class sorts the raw output of RefCov to something that is easier to interpret.
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('refcov_sort', '0.3')

    def _execute_tool(self):
        """
        Runs RefCovSort
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV key is required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        super(RefCovSort, self)._check_input()
        if 'TSV' not in self._tool_inputs or len(self._tool_inputs['TSV']) > 1:
                raise InvalidToolInputError('Invalid input given for RefCovSort (TSV and 1 file): {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Only TSV allowed as input for RefCovSort: {!r}'.format(self._tool_inputs))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs['TSV'][0].basename
        return os.path.join(self._folder, os.path.splitext(infile)[0])

    def __get_output_name(self):
        """
        Returns the name for the output file
        :return: Output file name
        """
        return self.__get_basename() + '_processed.tsv'

    def __set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(Path(self.__get_output_name()))]

    def __build_command(self):
        """
        Builds the sort and filter commands to process RefCov output
        :return: Command to run
        """
        if 'post_sort_by' in self._parameters:
            sort_by = self._parameters['post_sort_by'].value
        else:
            sort_by = 'coverage'
        if 'post_cov_cutoff' in self._parameters and 'post_depth_cutoff' in self._parameters:
            awk_cmd = "| awk '{{if($2>={0} && $6>={1})print;}}'".format(self._parameters['post_cov_cutoff'].value,
                                                                        self._parameters['post_depth_cutoff'].value)
        elif 'post_cov_cutoff' in self._parameters:
            awk_cmd = "| awk '{{if($2>={0})print;}}'".format(self._parameters['post_cov_cutoff'].value)
        elif 'post_depth_cutoff' in self._parameters:
            awk_cmd = "| awk '{{if($6>={0})print;}}'".format(self._parameters['post_depth_cutoff'].value)
        else:
            awk_cmd = ''
        # Sort by coverage (default, column 2) then sort by read depth (column 6)
        if sort_by == 'coverage':
            sort_cmd = "sort -k2,2nr -k6,6nr -t$'\t'"
        # Else: sort by read depth first (column 6) then sort by coverage (column 2)
        else:
            sort_cmd = "sort -k6,6nr -k2,2nr -t$'\t'"
        self._command.command = ' '.join([
            sort_cmd, str(self._tool_inputs['TSV'][0].path), awk_cmd, '>', self.__get_output_name()])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
