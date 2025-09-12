import re
import abc

from camel.app.command.command import Command
from camel.app.error import InvalidToolInputError, ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class GATK(Tool, metaclass=abc.ABCMeta):
    """
    Super class for GATK tools
    """

    def __init__(self, tool_name, version):
        """
        Initialize a GATK tool
        :param tool_name: Tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)
        self._specific_parameters = []
        self._required_inputs = []
        self._input_string = ''
        self._output_string = ''
        self._option_string = ''
        self._output_type = ''

    def _execute_tool(self):
        """
        Run a GATK function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidToolInputError('GATK {!r} required {!r} input is missing in _tool_inputs!'.format(
                    self._name, input_file))

        super(GATK, self)._check_input()

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

        if 'VCF' in self._tool_inputs:
            self._input_string += "-V {} ".format(self._tool_inputs['VCF'][0].path)

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs[self._output_type] = [ToolIOFile(self._folder / self.get_param_value('output'))]

    def _set_specific_parameters(self):
        """
        Set specific parameters that need special handling when required
        :return: None
        """
        pass

    def _build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._option_string += " ".join(self._build_options(excluded_parameters=self._specific_parameters))
        self._command.command = " ".join([
            self._tool_command, self._input_string, self._output_string, self._option_string
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Check the result of GATK tool run
        :param command: Command to check
        :return: None
        """
        if command.stdout == "":
            raise ToolExecutionError(self.name, "Fails to run as stdout is empty.")
        if not re.match('Exit status: 0', command.stdout.split('\n')[-2].rstrip()):
            raise ToolExecutionError(self.name, command.stdout)

    def _set_informs(self):
        """
        Set informs by analyze the output
        :return: None
        """
        # log WARNINGS in info.log
        for l in self._command.stdout.split('\n'):
            if re.match(r'WARNING', l):
                logger.info(" GATK - {}".format(l))

            # E.g., The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428
            match = re.search('The Genome Analysis Toolkit (GATK) (.+),', l)
            if match is not None:
                self.informs['tool_name'] = 'GATK'
                self.informs['tool_version'] = match.group(1)
            # Total filtering statistics
            #   0 reads were filtered out during the traversal out of approximately 1090654 total reads (0.00%)
            match = re.search(
                r'(\d+) reads were filtered out during the traversal out of approximately (\d+) total reads \((.+%)\)', l)
            if match is not None:
                self.informs['reads_total'] = match.group(2)
                if match.group(1) != '0':
                    self.informs['filtered_reads'] = "{}({})".format(match.group(1), match.group(3))
            # per Filter statistics:
            #    0 reads (0.00% of total) failing BadCigarFilter
            match = re.search(r'(\d+) reads \((.+%) of total\) failing (.+)', l)
            if match is not None and match.group(1) != '0':
                self.informs[match.group(3)] = "{}({})".format(match.group(1), match.group(2))
