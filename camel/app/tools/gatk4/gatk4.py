import abc
import re
import shutil
import tempfile

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class GATK4(Tool, metaclass=abc.ABCMeta):
    """
    Super class for GATK4 tools
    """

    def __init__(self, tool_name: str, version: str) -> None:
        """
        Initialize a GATK4 tool
        :param tool_name: Tool name
        :param version: Tool version
        :return: None
        """
        super().__init__(tool_name, version)
        self._specific_parameters = []
        self._required_inputs = []
        self._input_string = ''
        self._option_string = ''
        self._output_type = ''
        self._temp_dir = tempfile.mkdtemp(prefix='gatk4_', dir=config.dir_temp)
        self._java_options = f'"-mx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Djava.io.tmpdir={self._temp_dir}"'

    def update_java_options(self, java_options: str) -> None:
        """
        Returns the formatted java options of this tool.
        :return: Name
        """
        logger.info(f"Java options updated: '{java_options}'")
        self._java_options = f'"{java_options} -Djava.io.tmpdir={self._temp_dir}"'

    def _execute_tool(self) -> None:
        """
        Run a GATK function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        try:
            self._execute_command()
        finally:
            shutil.rmtree(self._temp_dir)
        self._set_informs()

    def _check_input(self) -> None:
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidToolInputError(f'GATK {self._name} required {input_file} input is missing in _tool_inputs!')
        super()._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += f"--intervals {self._tool_inputs['TXT_intervals'][0].path} "

        if 'VCF' in self._tool_inputs:
            self._input_string += f"--variant {self._tool_inputs['VCF'][0].path} "

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f"--reference {self._tool_inputs['FASTA_REF'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs[self._output_type] = [
            ToolIOFile(self.folder / self._parameters['output'].value)
        ]

    def _set_specific_parameters(self) -> None:
        """
        Set specific parameters that need special handling when required
        :return: None
        """
        pass

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        self._option_string += " ".join(self._build_options(excluded_parameters=self._specific_parameters))
        self._command.command = " ".join([
            'gatk --java-options', self._java_options, self._tool_command, self._input_string, self._option_string
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Check the result of the GATK run.
        :param command: Command to check
        :return: None
        """
        if self._command.exit_code != 0:
            raise ToolExecutionError(self.name, f'Error message:\n{command.stderr}')

    def _set_informs(self) -> None:
        """
        Set informs by analyze the output
        :return: None
        """
        # log WARNINGS in info.log
        for line in self._command.stdout.split('\n'):
            if re.match('WARN', line):
                logger.info(f" GATK - {line}")

            # E.g., The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428
            match = re.search('The Genome Analysis Toolkit (GATK) (.+),', line)
            if match is not None:
                self.informs['tool_name'] = 'GATK'
                self.informs['tool_version'] = match.group(1)
            # Total filtering statistics
            #   0 reads were filtered out during the traversal out of approximately 1090654 total reads (0.00%)
            match = re.search(
                r'(\d+) reads were filtered out during the traversal out of approximately (\d+) total reads \((.+%)\)', line)
            if match is not None:
                self.informs['reads_total'] = match.group(2)
                if match.group(1) != '0':
                    self.informs['filtered_reads'] = f"{match.group(1)}({match.group(3)})"
            # per Filter statistics:
            #    0 reads (0.00% of total) failing BadCigarFilter
            match = re.search(r'(\d+) reads \((.+%) of total\) failing (.+)', line)
            if match is not None and match.group(1) != '0':
                self.informs[match.group(3)] = f"{match.group(1)}({match.group(2)})"
