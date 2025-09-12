from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurClusterSplit(Mothur):
    """
    The cluster.split command can be used to assign sequences to OTUs and outputs a .list, .rabund, .sabund files.
    It splits large distance matrices into smaller pieces.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_cluster_split', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only the following keys are allowed: 'PHY', 'DIST', 'TSV_Names', 'TSV_Taxonomy', 'FASTA',
          'TSV_File', 'TSV_Counts'
        - Only one input file per key allowed
        - DIST, PHY and FASTA are mutually exclusive but this is not checked
        :return: None
        """
        super()._check_input()
        for key, input_files in self._tool_inputs.items():
            if key not in ['PHY', 'DIST', 'TSV_Names', 'TSV_Taxonomy', 'FASTA', 'TSV_File', 'TSV_Counts']:
                raise InvalidToolInputError('Invalid input key given for Mothur '
                                                     'cluster.split: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     cluster.split: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        # DIST, PHY and FASTA are mutually exclusive
        if 'DIST' in self._tool_inputs:
            items.append('column={}'.format(self._tool_inputs['fDIST'][0]))
        elif 'PHY' in self._tool_inputs:
            items.append('phylip={}'.format(self._tool_inputs['PHY'][0]))
        elif 'FASTA' in self._tool_inputs:
            items.append('fasta={}'.format(self._tool_inputs['FASTA'][0]))
        elif 'TSV_File' in self._tool_inputs:
            items.append('file={}'.format(self._tool_inputs['TSV_File'][0]))
        # Other keys are not mutually exclusive
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        if 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        if 'TSV_Taxonomy' in self._tool_inputs:
            items.append('taxonomy={}'.format(self._tool_inputs['TSV_Taxonomy'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = ''
        # DIST, PHY and FASTA are mutually exclusive
        if 'DIST' in self._tool_inputs:
            basename = super()._get_basename('DIST')
        elif 'PHY' in self._tool_inputs:
            basename = super()._get_basename('PHY')
        elif 'FASTA' in self._tool_inputs:
            basename = super()._get_basename()
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(basename + self.__get_extension())]

    def __get_extension(self):
        """
        Checks whether a different method is specified in the options as the ouput file names are based on the method
        specified. Also checks whether the cluster option is set to false as this gives another type of output.
        :return: String with extension
        """
        extension = ''
        for name, parameter in self._parameters.items():
            if name == 'method':
                if parameter.value == 'nearest':
                    extension = '.nn.unique_list.list'
                elif parameter.value == 'furthest':
                    extension = '.fn.unique_list.list'
            if name == 'cluster':
                if parameter.value.upper() == 'FALSE':
                    return '.file'
        return '.an.unique_list.list' if extension == '' else extension

    def __get_output_key(self):
        """
        Checks whether the cluster option is set to true as this
        changes the output.
        :return: Key to be used for the output
        """
        for name, parameter in self._parameters.items():
            if name == 'cluster':
                if parameter.value.upper() == 'FALSE':
                    return 'TSV_File'
                else:
                    return 'TSV_List'
        return 'TSV_List'
