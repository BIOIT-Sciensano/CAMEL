import abc
import os
import random
import tempfile
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.tool import Tool
from camel.app.loggers import logger


class Mothur(Tool):
    """
    A super class that contains all functions that are shared between
    Mothur commands.
    """

    def __init__(self, name, version=None):
        """
        Initialize tool
        :param name: Name of the tool
        :param version: Version of the tool
        :return: None
        """
        super().__init__(name, version)
        # For reproducibility a seed is specified for each operation
        self._seed = random.randint(1, 10000000)
        logger.debug(f'Set seed to: {self._seed}')
        self.__symlinks = []
        self.__temp_dir = tempfile.mkdtemp(dir='/scratch/temp')

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command.split()[0]} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split('=')[1].splitlines()[0].strip()

    def _execute_tool(self):
        """
        Runs Mothur
        :return: None
        """
        self._create_symlinks()
        self._build_command()
        self._execute_command()
        self._set_output()
        self._symlink_cleanup()

    def _create_symlinks(self):
        """
        Creates symlinks to all input files as Mothur does not allow the '-' character in the file names
        :return: None
        """
        new_inputs = {}
        for input_key, input_list in self._tool_inputs.items():
            new_inputs[input_key] = []
            for tool_input in input_list:
                if '-' in tool_input.path:
                    link_name = os.path.join(self.__temp_dir, tool_input.basename.replace('-', '_'))
                    os.symlink(tool_input.path, link_name)
                    self.__symlinks.append(link_name)
                    new_inputs[input_key] += [ToolIOFile(Path(link_name))]
                else:
                    new_inputs[input_key] += [tool_input]
        self._tool_inputs = new_inputs

    def _symlink_cleanup(self):
        """
        Removes all symlinks and the temporary directory that were created by the tool
        :return: None
        """
        for link in self.__symlinks:
            os.remove(link)
        os.rmdir(self.__temp_dir)

    def _build_command(self):
        """
        Concatenates required parameters and options to build the command to run
        :return: None
        """
        if self._tool_command.count('{') == 2:
            self._command.command = self._tool_command.format(self._build_input_string(), self._build_options())
        else:
            self._command.command = self._tool_command.format(self._build_input_string())

    @abc.abstractmethod
    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        :return: None
        """
        pass

    @abc.abstractmethod
    def _set_output(self):
        """
        Sets the name of the output files in the output file object
        :return: None
        """
        pass

    def _build_options(self, excluded_parameters=None, separator='='):
        """
        Creates the string with all the specified parameters
        :param excluded_parameters: list of parameters to be skipped (Optional)
        :param separator: separator used to combine the option and value (Optional)
        :return: String with command parameters
        """
        option_list = super()._build_options(excluded_parameters, separator)
        option_list += ['seed=' + str(self._seed)]
        return ', ' + ', '.join(option_list)

    def _get_labels(self):
        """
        Returns the labels that are in a list file or the ones that are specified as a parameter.
        :return: List of labels
        """
        if 'label' in self._parameters:
            return self._parameters['label'][1].strip().split('-')
        # If no label parameter is specified all the labels in the file will be used
        with open(self._tool_inputs['TSV_List'][0].path) as label_file:
            label_file.readline()
            return [line.split(None, 1)[0] for line in label_file]

    def _check_command_output(self, command: Command) -> None:
        """
        Analyzes output to discover if the run was successful. If an error was present in stdout, a RuntimeError is
        raised and stdout is displayed
        :return: None
        """
        for line in command.stdout.splitlines():
            if line.startswith('[ERROR]') or line.startswith('Unable to open'):
                raise RuntimeError('\n'.join(
                    command.stdout.splitlines()) + '\n' + '!!! Mothur failed to run !!! See above for more information.')

    def _get_basename(self, input_key='FASTA', suffix='.'):
        """
        Returns the prefix that will be used in the output. Example: Input file /test/data/file1.run1.fastq will return
        the following prefix: /test/data/file1.run1 (suffix = '.fastq')
        :param suffix: String that indicates the point where the path has to be cut (from the right)
        :param input_key: Key of the input file to be used
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs[input_key][0].basename
        return os.path.join(self._folder, infile[:infile.rfind(suffix)])

    def _get_extension(self, input_key='FASTA'):
        """
        Returns the extension of the file
        :param input_key: Key of the input file to be used
        :return: Extension of the file
        """
        return os.path.splitext(self._tool_inputs[input_key][0].path)[1]
