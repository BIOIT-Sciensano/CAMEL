from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurCluster(Mothur):
    """
    The dist.seqs command will calculate uncorrected pairwise distances between aligned DNA sequences. This approach
    is better than the commonly used DNADIST because the distances are not stored in RAM, rather they are printed
    directly to a file. Furthermore, it is possible to ignore "large" distances that one might not be interested in.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_cluster')
        self._required_input = []
        self._optional_input = ['TSV_Names', 'TSV_Counts']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid cfr. super.
        Additionally:
        - Either DIST or PHY or FASTA is required (not a combination)
        - Possible additional keys: TSV_Names, TSV_Counts (mutually exclusive)
        - When running this tool with the Vsearch method, FASTA and TSV_Counts are both required
        :return: None
        """
        # Check if one and only one of the required inputs is present
        if len([x for x in self._tool_inputs if x in ['FASTA', 'DIST', 'PHY']]) != 1:
            raise InvalidToolInputError("Invalid input files (keys) given. 'DIST', 'PHY' and 'FASTA' are mutually exclusive")
        # Add the present required input to the list
        self._required_input = [x for x in self._tool_inputs if x in ['FASTA', 'DIST', 'PHY']]
        # Check the mutually exclusive optional inputs
        if 'TSV_Names' in self._tool_inputs and 'TSV_Counts' in self._tool_inputs:
            raise InvalidToolInputError("Invalid input files (keys) given. 'TSV_Names' and 'TSV_Counts' are mutually exclusive")
        # Running with Vsearch method
        if self.get_param_value('method') in ['acg', 'dcg']:
            self._required_input = ['FASTA', 'TSV_Counts']
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        # Only DIST or PHY or FASTA is allowed, not more than one
        if 'DIST' in self._tool_inputs:
            items.append(f"column={self._tool_inputs['DIST'][0]}")
        elif 'PHY' in self._tool_inputs:
            items.append(f"phylip={self._tool_inputs['PHY'][0]}")
        elif 'FASTA' in self._tool_inputs:
            items.append(f"fasta={self._tool_inputs['FASTA'][0]}")
        # Either TSV_Counts or TSV_Names can be given, not both
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        elif 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        items.append(f"outputdir={self._folder}")
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # Only DIST or PHY is allowed, not both
        if 'DIST' in self._tool_inputs:
            basename = self._get_basename('DIST')
        elif 'PHY' in self._tool_inputs:
            basename = self._get_basename('PHY')
        else:
            basename = self._get_basename()
        self._tool_outputs['TSV_List'] = [ToolIOFile(basename.with_suffix(self.__get_extension()))]

    def __get_extension(self) -> str:
        """
        Checks whether a different method is specified in the options as
        the ouput file names are based on the method specified.
        :return: String with extension
        """
        method = self.get_param_value('method')
        if method == 'nearest':
            return '.nn.unique_list.list'
        elif method == 'furthest':
            return '.fn.unique_list.list'
        elif method == 'average':
            return '.an.unique_list.list'
        elif method == 'opti':
            return '.opti_mcc.list'
        elif method == 'agc':
            return '.agc.unique_list.list'
        elif method == 'dgc':
            return '.dgc.unique_list.list'
        raise ValueError(f'Invalid method: {method}')
