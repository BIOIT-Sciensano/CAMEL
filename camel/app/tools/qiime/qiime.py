import abc
import os

from camel.app.command.command import Command
from camel.app.error import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class Qiime(Tool):
    """
    A super class that contains all functions that are shared between Qiime scripts.
    """

    def __init__(self, name, version):
        """
        Initialize tool
        :param name: Name of the tool
        :param version: Version of the tool
        :return: None
        """
        super().__init__(name, version)
        self._parameter_file = 'parameters.txt'

    def _execute_tool(self):
        """
        Runs Prinseq
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
        self._add_log_filename()

    def _build_options(self, excluded_parameters=None, separator=' '):
        """
        Creates a file with parameters that are not part of the parameter string and returns the parameter string to
        be used.
        :param excluded_parameters: list of parameters to be skipped (Optional)
        :param separator: separator used to combine the option and value (Optional)
        :return: String with command parameters
        """
        file_params = []
        for name, param in self._parameters.items():
            if ':' in name:
                self.__write_to_parameter_file(param)
                file_params.append(name)
        if excluded_parameters is not None:
            file_params += excluded_parameters
        return super(Qiime, self)._build_options(excluded_parameters=file_params, delimiter=separator)

    def _build_command(self):
        """
        Concatenates required parameters and options to build the command to run
        :return: None
        """
        # _build_options needs to be executed first as _build_input_string will use the parameter file
        options = self._build_options()
        self._command.command = '{} {} {}'.format(self._tool_command, self._build_input_string(), options)

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

    def __write_to_parameter_file(self, parameter):
        """
        Write the given parameter to the parameter file
        :param parameter: Paramater dictionary to write to the file
        :return: None
        """
        param_file = self._folder / self._parameter_file
        with open(param_file, 'a') as outfile:
            outfile.write(parameter.option + ' ' + parameter.value + '\n')

    def _get_basename(self, input_key='FASTA', suffix='.'):
        """
        Returns the prefix that will be used in the output.
        Example: Input file /test/data/file1.run1.fastq will return the
        following prefix: /test/data/file1.run1 (suffix = '.fastq')
        :param suffix: Suffix that has to be removed
        :param input_key: Key of the input file to be used
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs[input_key][0].path)
        return os.path.join(self._folder, infile[:infile.rfind(suffix)])

    def _add_log_filename(self):
        """
        Searches for the log file that was created during the run and adds this to the output
        :return: None
        """
        log_files = [ToolIOFile(self._folder / f) for f in os.listdir(self._folder) if f.startswith('log_')]
        if len(log_files) > 0:
            self._tool_outputs['LOG'] = log_files
            logger.debug('Added log file to outputs: {}'.format(self._tool_outputs))

    def _check_command_output(self, command: Command):
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if command.stderr != '' and 'error' in command.stderr.lower():
            raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        elif command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Command execution failed (Exit code: {self._command.exit_code})")
