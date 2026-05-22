from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurGetGroups(Mothur):
    """
    The get.groups command selects sequences from a specific group or set of groups from the following file types:
    fasta, name, group, list, taxonomy.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_get_groups', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Either TSV_Groups or TSV_Counts are required
        - Additional allowed keys are: 'FASTA', 'TSV_Names', 'TSV_Accnos', 'TSV_List', 'TSV_Taxonomy'
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        if 'TSV_Groups' not in self._tool_inputs and 'TSV_Counts' not in self._tool_inputs:
            raise InvalidToolInputError('Not enough valid input files given for Mothur '
                                                 f'get.groups: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
                           'TSV_Accnos', 'TSV_List', 'TSV_Taxonomy']:
                raise InvalidToolInputError(f'Invalid input key given for Mothur get.groups: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files given for Mothur \
                                                     get.groups: {self._tool_inputs!r}')

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        input_parameters = {'FASTA': 'fasta=',
                            'TSV_Names': 'name=',
                            'TSV_Counts': 'count=',
                            'TSV_Groups': 'group=',
                            'TSV_List': 'list=',
                            'TSV_Accnos': 'accnos=',
                            'TSV_Taxonomy': 'taxonomy='}
        # Based on the key the correct option flag is added to the input string
        for key, input_files in self._tool_inputs.items():
            items.append(f'{input_parameters[key]}{input_files[0].path}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        output_extensions = {'FASTA': ['.', '.pick.fasta'],
                             'TSV_Names': ['.', '.pick.names'],
                             'TSV_Counts': ['.', '.pick.count_table'],
                             'TSV_Groups': ['.', '.pick.groups'],
                             'TSV_List': ['.', '.pick.list'],
                             'TSV_Taxonomy': ['.', '.pick.taxonomy']}
        # Based on the key the correct output file is added to the input string
        for key, input_files in self._tool_inputs.items():
            basename = super()._get_basename(key, output_extensions[key][0])
            self._tool_outputs[key] = [ToolIOFile(basename + output_extensions[key][1])]
