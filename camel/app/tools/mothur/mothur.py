import abc
import random
import tempfile
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class Mothur(Tool, metaclass=abc.ABCMeta):
    """
    A super class that contains all functions that are shared between Mothur commands.
    """

    def __init__(self, tool_name: str) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__(tool_name, version=None)
        # For reproducibility a seed is specified for each operation
        self._seed = random.randint(1, 10000000)
        logger.info(f'Set seed to: {self._seed}')
        self._temp_dir = Path(tempfile.mkdtemp(prefix='mothur_', dir=config.dir_temp))
        self._required_input = []
        self._optional_input = []

    def get_version(self) -> str:
        """
        Retrieves the mothur version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command.split()[0]} --version')
        self._execute_command(command, is_version_cmd=True)
        for line in command.stdout.splitlines():
            if 'Mothur version=' in line:
                return line.split('=')[1].strip()
        raise RuntimeError(f'Failed to parse mothur version from: {command.stdout}')

    def _execute_tool(self) -> None:
        """
        Runs Mothur.
        :return: None
        """
        self._create_symlinks(self._temp_dir)
        self._build_command()
        self._execute_command()
        self._set_output()

    def _create_symlinks(self, dir_temp: Path) -> None:
        """
        Creates symlinks to all input files as Mothur does not allow the '-' character in the file names.
        :return: None
        """
        new_inputs = {}
        for input_key, input_list in self._tool_inputs.items():
            new_inputs[input_key] = []
            for tool_input in input_list:
                if '-' in str(tool_input.path):
                    link_name = dir_temp / tool_input.basename.replace('-', '_')
                    link_name.symlink_to(tool_input.path)
                    new_inputs[input_key] += [ToolIOFile(link_name)]
                else:
                    new_inputs[input_key] += [tool_input]
        self._tool_inputs = new_inputs

    def _build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command to run
        :return: None
        """
        if self._tool_command.count('{') == 2:
            self._command.command = self._tool_command.format(self._build_input_string(), self._build_options())
        else:
            self._command.command = self._tool_command.format(self._build_input_string())

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        # Check if all required inputs are present
        for key in self._required_input:
            if key not in self._tool_inputs:
                raise InvalidToolInputError(f"Input '{key}' is required")
        # Check whether all present input files are correct and there is only one file per key
        for tool_input in self._tool_inputs:
            if tool_input not in self._required_input + self._optional_input:
                raise InvalidToolInputError(f'Invalid input key: {tool_input}')
            if tool_input in ['FASTQ_PE', 'FASTA_PE']:
                if len(self._tool_inputs[tool_input]) != 2:
                    raise InvalidToolInputError(f'Invalid number of files given for Mothur make.contigs: {tool_input}')
            elif len(self._tool_inputs[tool_input]) != 1:
                raise InvalidToolInputError(f"Invalid number (max = 1) of files per key given for '{tool_input}'.")
        super()._check_input()

    @abc.abstractmethod
    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories.
        :return: Input string
        """
        pass

    @abc.abstractmethod
    def _set_output(self) -> None:
        """
        Sets the name of the output files in the output file object.
        :return: None
        """
        pass

    def _build_options(self, excluded_parameters: list[str] | None = None, separator: str = '=') -> str:
        """
        Creates the string with all the specified parameters.
        :param excluded_parameters: list of parameters to be skipped (Optional)
        :param separator: separator used to combine the option and value (Optional)
        :return: String with command parameters
        """
        option_list = super()._build_options(excluded_parameters, separator)
        option_list.append(f'seed={self._seed}')
        return ', ' + ', '.join(option_list)

    def _get_labels(self) -> list[str]:
        """
        Returns the labels that are in a list file or the ones that are specified as a parameter.
        :return: List of labels
        """
        if 'label' in self._parameters:
            return self.get_param_value('label').strip().split('-')
        # If no label parameter is specified all the labels in the file will be used
        with open(self._tool_inputs['TSV_List'][0].path) as label_file:
            label_file.readline()
            return [line.split(None, 1)[0] for line in label_file]

    def _get_basename(self, input_key: str = 'FASTA', suffixes_to_remove: set[str] = None) -> Path:
        """
        Returns the basename that will be used in the output. Example: Input file /test/data/file1.run1.fastq will return
        the following prefix: /output/directory/file1.run1.
        :param input_key: Key of the input file to be used
        :param suffixes_to_remove: Set of suffixes that need to be removed from the path
        :return: String with the prefix used in the output
        """
        basename = self._folder / self._tool_inputs[input_key][0].basename
        if suffixes_to_remove is None:
            return basename
        else:
            while basename.suffix in suffixes_to_remove:
                basename = basename.with_suffix('')
            return basename

    def _get_extension(self, input_key: str = 'FASTA') -> str:
        """
        Returns the extension of the file
        :param input_key: Key of the input file to be used
        :return: Extension of the file
        """
        return self._tool_inputs[input_key][0].path.suffix
