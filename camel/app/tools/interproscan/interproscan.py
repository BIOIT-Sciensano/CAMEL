import os
import shutil
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class Interproscan(Tool):
    """
    InterPro is a database which integrates together predictive information about proteins' function from a number of
    partner resources, giving an overview of the families that a protein belongs to and the domains and sites it
    contains. InterProScan is the software that can be used to query this database
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('interproscan', '5.20-59.0')
        self.__input_key = None

    def _execute_tool(self):
        """
        Runs InterProScan
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        shutil.rmtree(self.__get_temp_dir())

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA_Nucl or FASTA_Prot is required
        - Only one input file allowed per input key
        :return: None
        """
        super(Interproscan, self)._check_input()
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Invalid number of input keys given for InterProScan, only FASTA_Nucl '
                                                 'or FASTA_Prot allowed: {!r}'.format(self._tool_inputs))
        if 'FASTA_Nucl' not in self._tool_inputs and 'FASTA_Prot' not in self._tool_inputs:
            raise InvalidToolInputError('No valid input key given for InterProScan, FASTA_Nucl or FASTA_Prot '
                                                 'needed: {!r}'.format(self._tool_inputs))
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidToolInputError('Too many input files per key given for InterProScan '
                                                     '(max = 1): {!r}'.format(self._tool_inputs))

    def __set_input_key(self):
        """
        Sets the input key that is present in the inputs
        :return: None
        """
        self.__input_key = 'FASTA_Nucl' if 'FASTA_Nucl' in self._tool_inputs else 'FASTA_Prot'

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        if 'formats' in self._parameters:
            output_keys = self._parameters['formats'][1].split(',')
        else:
            output_keys = ['GFF3', 'XML', 'TSV']
        for key in output_keys:
            self._tool_outputs[key] = [ToolIOFile(Path(
                f'{os.path.join(self._folder, self.__get_basename())}.{key.lower()}'))]

    def __get_basename(self):
        """
        Creates the basename for the output files
        :return: Basename
        """
        return self._tool_inputs[self.__input_key][0].basename

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        seqtype = 'n' if self.__input_key == 'FASTA_Nucl' else 'p'
        return '--input {} --seqtype {}'.format(self._tool_inputs[self.__input_key][0], seqtype)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = self.__build_options_string()
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __get_temp_dir(self):
        """
        Returns the path to the temporary directory and creates it if it does not yet exist
        :return: Path to temporary directory
        """
        if os.path.isfile(os.path.join(self._folder, 'temp_ips')):
            raise IOError('A file with the name temp_ips already exists!')
        if not os.path.isdir(os.path.join(self._folder, 'temp_ips')):
            os.mkdir(os.path.join(self._folder, 'temp_ips'))
        return os.path.join(self._folder, 'temp_ips')

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        if 'Unfortunately the web service has failed.' in command.stdout:
            logger.warning('The local lookup webservice was not reachable by InterProScan!')
        if command.exit_code != 0:
            raise ToolExecutionError(
                self.name, f"Command execution failed for InterProScan (Exit code: {command.exit_code})")

    def __build_options_string(self):
        """
        Creates the string with all the specified parameters
        :return: String with command parameters
        """
        option_list = super(Interproscan, self)._build_options()
        option_list.append('--output-dir {}'.format(self._folder))
        option_list.append('--tempdir {}'.format(self.__get_temp_dir()))
        return ' '.join(option_list)
