from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurCluster(Mothur):
    """
    The dist.seqs command will calculate uncorrected pairwise distances between aligned DNA sequences. This approach
    is better than the commonly used DNADIST because the distances are not stored in RAM, rather they are printed
    directly to a file. Furthermore, it is possible to ignore "large" distances that one might not be interested in.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_cluster', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Either DIST or PHY or FASTA is required (not a combination)
        - Possible additional keys: TSV_Names, TSV_Counts (mutually exclusive)
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        if len([x for x in self._tool_inputs if x in ['FASTA', 'DIST', 'PHY']]) != 1:
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur cluster, only DIST or '
                                                 'PHY or FASTA allowed: {!r}'.format(self._tool_inputs))
        if 'TSV_Names' in self._tool_inputs and 'TSV_Counts' in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur cluster, TSV_Names and '
                                                 'TSV_Counts not allowed together: {!r}'.format(self._tool_inputs))
        if self._parameters['method'].value in ['acg', 'dcg']:
            if 'FASTA' not in self._tool_inputs or 'TSV_Counts' not in self._tool_inputs:
                raise InvalidToolInputError('Fasta and count table required for Mothur cluster with Vsearch '
                                                     'method: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['DIST', 'TSV_Names', 'TSV_Counts', 'PHY', 'FASTA']:
                raise InvalidToolInputError('Invalid input key given for Mothur '
                                                     'cluster: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     cluster: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        # Only DIST or PHY is allowed, not both
        if 'DIST' in self._tool_inputs:
            items.append('column={}'.format(self._tool_inputs['DIST'][0]))
        elif 'PHY' in self._tool_inputs:
            items.append('phylip={}'.format(self._tool_inputs['PHY'][0]))
        elif 'FASTA' in self._tool_inputs:
            items.append('fasta={}'.format(self._tool_inputs['FASTA'][0]))
        # Either TSV_Counts or TSV_Names can be given, not both
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
        basename = ''
        # Only DIST or PHY is allowed, not both
        if 'DIST' in self._tool_inputs:
            basename = super()._get_basename('DIST')
        elif 'PHY' in self._tool_inputs:
            basename = super()._get_basename('PHY')
        elif 'FASTA' in self._tool_inputs:
            basename = super()._get_basename()
        self._tool_outputs['TSV_List'] = [ToolIOFile(basename + self.__get_extension())]

    def __get_extension(self):
        """
        Checks whether a different method is specified in the options as
        the ouput file names are based on the method specified.
        :return: String with extension
        """
        parameter = self._parameters['method']
        if parameter.value == 'nearest':
            return '.nn.unique_list.list'
        elif parameter.value == 'furthest':
            return '.fn.unique_list.list'
        elif parameter.value == 'average':
            return '.an.unique_list.list'
        elif parameter.value == 'opti':
            return '.opti_mcc.list'
        elif parameter.value == 'agc':
            return '.agc.unique_list.list'
        elif parameter.value == 'dgc':
            return '.dgc.unique_list.list'
