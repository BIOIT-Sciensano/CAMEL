import errno
import shutil

import os
import re

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Phymmbl(Tool):
    """
    Phymm, a new classification approach for metagenomics data which uses interpolated Markov models (IMMs) to
    taxonomically classify DNA sequences, can accurately classify reads as short as 100 bp. Its accuracy for short
    reads represents a significant leap forward over previous composition-based classification methods. PhymmBL, the
    hybrid classifier included in this distribution which combines analysis from both Phymm and BLAST, produces even
    higher accuracy.
    """

    def __init__(self):
        """
        Initialize tool
                :return: None
        """
        super().__init__('phymmbl', '4.0')

    def _execute_tool(self):
        """
        Function to run tool
        :return: None
        """
        self.__build_command()
        self.__create_phymmbl_copy()
        self._execute_command(dir_=self._folder / 'phymmbl')
        self.__clean_temp_directory()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only FASTA and DB are allowed and required
        - Only one input file allowed
        :return: None
        """
        super(Phymmbl, self)._check_input()
        if 'FASTA' not in self._tool_inputs or 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input key given for PhymmBL, '
                                                 'only FASTA and DB allowed: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 2:
            raise InvalidToolInputError('Invalid number of input keys given for PhymmBL, '
                                                 'only FASTA and DB allowed: {!r}'.format(self._tool_inputs))
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files per '
                                                     'key given for PhymmBL: {!r}'.format(self._tool_inputs))

    def __build_output_name(self, prefix):
        """
        Builds the complete name of the output file
        :param prefix: Prefix to use
        :return: Complete path to the output file
        """
        basename = prefix + re.sub(r'[/.]', '_', str(self._tool_inputs['FASTA'][0].path)) + '.txt'
        return self._folder / 'phymmbl' / basename

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['TSV_RawPhymm'] = [ToolIOFile(self.__build_output_name('rawPhymmOutput_'))]
        self._tool_outputs['TSV_RawBLAST'] = [ToolIOFile(self.__build_output_name('rawBlastOutput_'))]
        self._tool_outputs['TSV_Phymm'] = [ToolIOFile(self.__build_output_name('results.01.phymm_'))]
        self._tool_outputs['TSV_BLAST'] = [ToolIOFile(self.__build_output_name('results.02.blast_'))]
        self._tool_outputs['TSV_PhymmBL'] = [ToolIOFile(self.__build_output_name('results.03.phymmBL_'))]

    def __create_phymmbl_copy(self):
        """
        Creates a local clone of the phymmbl binaries and links all necessary database directories + a copy of the input
        :return: None
        """
        shutil.copytree(self.__get_phymmbl_path(), os.path.join(self._folder, 'phymmbl'))
        self.__create_symlinks()

    def __get_phymmbl_path(self):
        """
        Returns the path to the phymmbl binaries by retrieving it from the system path after loading the relevant
        module.
        :return: Path to phymmbl binaries
        """
        phymm_cmd = Command()
        phymm_cmd.command = self._build_dependencies() + ' echo $PATH'
        phymm_cmd.run_command(self._folder)
        for item in phymm_cmd.stdout.split(':'):
            if '/phymmbl/' in item:
                return item
        raise OSError('The PhymmBL executable location was not found. Is the correct Lmod module loaded?')

    def __create_symlinks(self):
        """
        Creates symbolic links for all the directories that are associated with the PhymmBL database
        :return: None
        """
        for directory in ['.genomeData', '.blastData', '.taxonomyData']:
            try:
                os.symlink(os.path.join(self._tool_inputs['DB'][0].path, directory),
                           os.path.join(self._folder, 'phymmbl', directory))
            except OSError as e:
                if e.errno == errno.EEXIST:
                    shutil.rmtree(os.path.join(self._folder, 'phymmbl', directory))
                    os.symlink(os.path.join(self._tool_inputs['DB'][0].path, directory),
                               os.path.join(self._folder, 'phymmbl', directory))
                else:
                    raise e

    def __clean_temp_directory(self):
        """
        Removes all non-output files from the running directory (i.e. all files that were copied from the binaries
        directory and symlinks that were created)
        :return: None
        """
        output_folder = os.path.join(self._folder, 'phymmbl')
        for item in os.listdir(output_folder):
            if item.startswith(('rawBlastOutput', 'rawPhymmOutput', 'results')):
                continue
            elif os.path.islink(os.path.join(output_folder, item)):
                os.unlink(os.path.join(output_folder, item).rstrip('/'))
            elif os.path.isdir(os.path.join(output_folder, item)):
                shutil.rmtree(os.path.join(output_folder, item))
            else:
                os.remove(os.path.join(output_folder, item))

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self._tool_inputs['FASTA'][0].path
        self._command.command = self._tool_command + ' {}'.format(input_string)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
