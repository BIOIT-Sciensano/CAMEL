import os
import re
import shutil

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Orphelia(Tool):
    """
    Orphelia is a metagenomic ORF finding tool for the prediction of protein coding genes in short, environmental DNA
    sequences with unknown phylogenetic origin
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('orphelia', '74')

    def _execute_tool(self):
        """
        Runs Orphelia
        :return: None
        """
        self.__build_command()
        self.__create_orphelia_copy()
        self._command.command = self.__remove_lmod(self._build_dependencies()) + self._command.command
        self._command.run_command(self.__get_working_folder())
        self._check_command_output(self._command)
        self.__clean_temp_directory()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only FASTA is allowed and required
        - Only one input file allowed
        :return: None
        """
        super(Orphelia, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError(f'Invalid input key given for Orphelia, only FASTA allowed: {self._tool_inputs!r}')
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError(f'Invalid number of input keys given for Orphelia, only FASTA allowed: {self._tool_inputs!r}')
        if len(self._tool_inputs['FASTA']) > 1:
            raise ValueError(f'Invalid number (max = 1) of files per key given for Orphelia: {self._tool_inputs!r}')

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['COORD'] = [ToolIOFile(self._folder / 'gene.pred')]
        self._tool_outputs['TEXT_Header'] = [ToolIOFile(self._folder / 'frags.header')]
        self._tool_outputs['TEXT_Seq'] = [ToolIOFile(self._folder / 'orf.seq')]
        self._tool_outputs['TEXT_Coords'] = [ToolIOFile(self._folder / 'orf.coords')]
        self._tool_outputs['TEXT_Tis'] = [ToolIOFile(self._folder / 'tis.seq')]

    def __create_orphelia_copy(self):
        """
        Creates a local clone of the orphelia binaries
        :return: None
        """
        shutil.copytree(self.__get_orphelia_path(), self.__get_working_folder())

    def __get_orphelia_path(self):
        """
        Returns the path to the orphelia binaries by retrieving it from the system path after loading the relevant
        module.
        :return: Path to orphelia binaries
        """
        orphelia_cmd = Command()
        orphelia_cmd.command = self._build_dependencies() + ' echo $PATH'
        orphelia_cmd.run(self._folder)
        for item in orphelia_cmd.stdout.split(':'):
            if '/orphelia/' in item:
                return item
        raise OSError('The Orphelia executable location was not found. Is the correct Lmod module loaded?')

    def __clean_temp_directory(self):
        """
        Removes the temporary version of Orphelia.
        :return: None
        """
        shutil.rmtree(self.__get_working_folder())

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self._tool_inputs['FASTA'][0].path
        self._tool_command = os.path.join(self.__get_working_folder(), self._tool_command)
        self._command.command = '{} -s {} -o {} {}'.format(self._tool_command, input_string, self._folder,
                                                           ' '.join(self._build_options()))

    @staticmethod
    def __remove_lmod(command):
        """
        Remove the orphelia module from the lmod load command
        :param command: Lmod load command
        :return: Lmod load command without the orphelia module (or empty string if no other modules are left)
        """
        sub = re.sub(r'orphelia/[0-9.\s]+', '', command)
        # A module load command should be at least 15 characters to actually have a module that needs to be loaded
        return sub if len(sub) > 15 else ''

    def __get_working_folder(self):
        """
        Returns the folder where the Orphelia clone is located
        :return: Folder where Orphelia clone is located
        """
        return os.path.join(self._folder, 'orphelia')

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
