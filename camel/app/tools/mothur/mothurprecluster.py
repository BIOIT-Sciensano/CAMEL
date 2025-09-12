from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.mothur.mothur import Mothur


class MothurPreCluster(Mothur):
    """
    The pre.cluster command implements a pseudo-single linkage algorithm with the goal of removing sequences that are
    likely due to pyrosequencing errors.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_pre_cluster', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - TSV_Names or TSV_Counts is required
        - TSV_Groups is an additional allowed key
        - Only one input file per key is allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs or ('TSV_Names' in self._tool_inputs and 'TSV_Counts' in self._tool_inputs):
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur '
                                                 'pre.cluster: {!r}'.format(self._tool_inputs))
        if 'TSV_Names' not in self._tool_inputs and 'TSV_Counts' not in self._tool_inputs:
            raise InvalidToolInputError('Missing input files (key) for Mothur '
                                                 'pre.cluster: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Counts', 'TSV_Names', 'TSV_Groups']:
                raise InvalidToolInputError('Invalid input key given for Mothur pre.cluster: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     pre.cluster: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0])]
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        elif 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        if 'skip_step' in self._parameters:
            self._tool_outputs['FASTA'] = [ToolIOFile(self._tool_inputs['FASTA'][0].path)]
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(self._tool_inputs['TSV_Counts'][0].path)]
        else:
            basename = super()._get_basename()
            self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.precluster.fasta')]
            group_names = []
            if 'TSV_Counts' in self._tool_inputs:
                group_names = self.__get_group_names()
                # One file is always in the output when TSV_Counts is specified (.precluster.count_table)
                self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename + '.precluster.count_table')]
            elif 'TSV_Names' in self._tool_inputs:
                group_names = self.__get_group_names()
            if len(group_names) != 0:
                self._tool_outputs['TSV_Map'] = []
                for item in group_names:
                    # For each group a precluster map file is created
                    self._tool_outputs['TSV_Map'].append(ToolIOFile(basename + '.precluster.' + item + '.map'))

    def __get_group_names(self):
        """
        Returns the names of the different groups that are in a count file. A names file can also be passed but this
        has to be implemented as it is not clear what the structure of such a file is at this moment.
        :return: List of group names
        """
        group_names = []
        if 'TSV_Counts' in self._tool_inputs:
            with open(self._tool_inputs['TSV_Counts'][0].path, 'r') as count_file:
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

    def _build_options(self, excluded_parameters=None, separator='='):
        """
        Creates the string with all the specified parameters
        :param excluded_parameters: list of parameters to be skipped (Optional)
        :param separator: separator used to combine the option and value (Optional)
        :return: String with command parameters
        """
        return super()._build_options(excluded_parameters=['skip_step'], separator=separator)

    def _execute_tool(self):
        """
        Runs Mothur Pre.cluster
        :return: None
        """
        if 'skip_step' not in self._parameters:
            self._create_symlinks()
            self._build_command()
            self._execute_command()
            self._set_output()
            self._symlink_cleanup()
        else:
            self._set_output()
            logger.warning("Skipping the precluster step as requested by setting the skip_step parameter!")
