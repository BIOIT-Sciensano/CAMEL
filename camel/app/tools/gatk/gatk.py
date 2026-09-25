import abc
import re
import tempfile
from pathlib import Path

from camelcore.app.command import Command

from camel.app.config import config
from camel.app.core.errors import ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class GATK(Tool, metaclass=abc.ABCMeta):
    """
    Base class for GATK (v3) tools.
    """

    def __init__(self, tool_name: str) -> None:
        """
        Initializes a GATK tool.
        :param tool_name: Tool name
        :return: None
        """
        super().__init__(tool_name, version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        # The GATK wrapper script writes a log file to the working directory
        with tempfile.TemporaryDirectory(dir=config.dir_temp) as dir_temp:
            self._execute_command(command, dir_=Path(dir_temp), is_version_cmd=True)
        return command.stdout.splitlines()[0].strip()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
        self.__set_informs()

    def _build_command(self) -> None:
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._get_input_options(),
            *self._get_parameter_options()
        ])

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files shared by the GATK tools.
        Subclasses extend this with their tool-specific inputs.
        :return: Input options
        """
        options = []
        if 'FASTA_REF' in self._tool_inputs:
            options.append(f"-R {self._tool_inputs['FASTA_REF'][0].path}")
        for interval_file in self._tool_inputs.get('TXT_intervals', []):
            options.append(f'-L {interval_file.path}')
        return options

    def _get_parameter_options(self) -> list[str]:
        """
        Returns the command line options for the parameters.
        Subclasses override this for parameters that need special handling.
        :return: Parameter options
        """
        return self._build_options()

    @abc.abstractmethod
    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        The GATK wrapper script merges stderr into stdout and exits with a non-zero code when an error is logged.
        :param command: Command to check
        :return: None
        """
        if command.exit_code == 0:
            return
        lines = command.stdout.splitlines()
        error_lines = [line for line in lines if line.startswith('##### ERROR MESSAGE')]
        raise ToolExecutionError(self.name, '\n'.join(error_lines if len(error_lines) > 0 else lines[-10:]))

    def __set_informs(self) -> None:
        """
        Sets the informs by parsing the tool output.
        :return: None
        """
        for line in self._command.stdout.splitlines():
            if line.startswith('WARN'):
                logger.info(f'GATK - {line}')

            # E.g., 0 reads were filtered out during the traversal out of approximately 1090654 total reads (0.00%)
            match = re.search(
                r'(\d+) reads were filtered out during the traversal out of approximately (\d+) total reads \((.+%)\)',
                line)
            if match is not None:
                self.informs['reads_total'] = match.group(2)
                if match.group(1) != '0':
                    self.informs['filtered_reads'] = f'{match.group(1)}({match.group(3)})'

            # E.g., 0 reads (0.00% of total) failing BadCigarFilter
            match = re.search(r'(\d+) reads \((.+%) of total\) failing (.+)', line)
            if match is not None and match.group(1) != '0':
                self.informs[match.group(3)] = f'{match.group(1)}({match.group(2)})'
