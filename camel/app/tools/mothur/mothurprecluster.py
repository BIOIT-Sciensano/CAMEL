import re

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.loggers import logger
from camel.app.tools.mothur.mothur import Mothur


class MothurPreCluster(Mothur):
    """
    The pre.cluster command implements a pseudo-single linkage algorithm with the goal of removing sequences that are
    likely due to pyrosequencing errors.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_pre_cluster')
        self._required_input = ['FASTA']
        self._optional_input = ['TSV_Groups']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - TSV_Names or TSV_Counts is required
        - TSV_Groups is an additional allowed key
        - Only one input file per key is allowed
        :return: None
        """
        if 'TSV_Counts' in self._tool_inputs:
            self._required_input.append('TSV_Counts')
        elif 'TSV_Names' in self._tool_inputs:
            self._required_input.append('TSV_Names')
        else:
            raise InvalidToolInputError('Either TSV_Counts or TSV_Names is required')
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        elif 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        if 'skip_step' in self._parameters:
            self._tool_outputs['FASTA'] = [ToolIOFile(self._tool_inputs['FASTA'][0].path)]
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(self._tool_inputs['TSV_Counts'][0].path)]
        else:
            basename = self._get_basename()
            self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix('.precluster.fasta'))]
            group_names = []
            if 'TSV_Counts' in self._tool_inputs:
                group_names = self.__get_group_names()
                # One file is always in the output when TSV_Counts is specified (.precluster.count_table)
                self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename.with_suffix('.precluster.count_table'))]
            elif 'TSV_Names' in self._tool_inputs:
                group_names = self.__get_group_names()
            if len(group_names) != 0:
                self._tool_outputs['TSV_Map'] = []
                for item in group_names:
                    # For each group a precluster map file is created
                    self._tool_outputs['TSV_Map'].append(ToolIOFile(basename.with_suffix(f'.precluster.{item}.map')))

    def __get_group_names(self) -> list[str]:
        """
        Returns the names of the different groups that are in a count file. A names file can also be passed but this
        has to be implemented as it is not clear what the structure of such a file is at this moment.
        :return: List of group names
        """
        group_names = []
        if 'TSV_Counts' in self._tool_inputs:
            with open(self._tool_inputs['TSV_Counts'][0].path) as count_file:
                # The first line contains the groups
                group_line = count_file.readline().strip()
            groups = group_line.split('\t')
            # The group names start from the third element
            for group in groups[2:]:
                group_names.append(group)
        # TSV_Names is not yet implemented
        elif 'TSV_Names' in self._tool_inputs:
            raise RuntimeError('Using a names file is not yet implemented for pre.cluster!')
        return group_names

    def _build_options(self, excluded_parameters: list[str] | None = None, separator: str = '=') -> str:
        """
        Creates the string with all the specified parameters
        :param excluded_parameters: list of parameters to be skipped (Optional)
        :param separator: separator used to combine the option and value (Optional)
        :return: String with command parameters
        """
        return super()._build_options(excluded_parameters=['skip_step'], separator=separator)

    def _execute_tool(self) -> None:
        """
        Runs Mothur Pre.cluster
        :return: None
        """
        if 'skip_step' not in self._parameters:
            if not self._count_file_has_groups():
                self._add_dummy_group()
            self._create_symlinks(self._temp_dir)
            self._build_command()
            self._execute_command()
            self._set_output()
        else:
            self._set_output()
            logger.warning("Skipping the precluster step as requested by setting the skip_step parameter!")

    def _count_file_has_groups(self) -> bool:
        """
        Checks whether the given counts file contains groups (starting from third column) as pre.cluster fails if this is not the case
        :return: True/False
        """
        with self._tool_inputs['TSV_Counts'][0].path.open() as inhandle:
            line = inhandle.readline()
        return len(re.split('[\t ]', line)) > 2

    def _add_dummy_group(self) -> None:
        """
        Adds a dummy group column to the counts file
        :return: None
        """
        command = f"awk '{{print $0,$NF}}' {self._tool_inputs['TSV_Counts'][0].path} | sed \"s/total total/total A/\""
        cmd = Command(command)
        cmd.run(self._folder, disable_logging=True)
        new_counts = self._folder / self._tool_inputs['TSV_Counts'][0].basename
        with new_counts.open('w') as outhandle:
            outhandle.write(cmd.stdout)
        self._tool_inputs['TSV_Counts'] = [ToolIOFile(new_counts)]
