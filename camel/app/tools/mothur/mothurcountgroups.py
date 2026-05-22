from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurCountGroups(Mothur):
    """
    The count.groups command counts sequences from a specific group or set of groups.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_count_groups', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only the following keys are allowed: 'TSV_Groups', 'TSV_Counts', 'TSV_Shared', 'TSV_Accnos'
        - TSV_Groups, TSV_Shared and TSV_Counts are mutually exclusive but this is not checked
        - Only one input file per key is allowed
        :return: None
        """
        super()._check_input()
        for key, input_files in self._tool_inputs.items():
            if key not in ['TSV_Groups', 'TSV_Counts', 'TSV_Shared', 'TSV_Accnos']:
                raise InvalidToolInputError('Invalid input key given for Mothur '
                                                     'count.groups: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     count.groups: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        elif 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        elif 'TSV_Shared' in self._tool_inputs:
            items.append('shared={}'.format(self._tool_inputs['TSV_Shared'][0]))
        if 'TSV_Accnos' in self._tool_inputs:
            items.append('reftaxonomy={}'.format(self._tool_inputs['TSV_Accnos'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them.
        :return: None
        """
        basename = ''
        # TSV_Groups, TSV_Shared and TSV_Counts are mutually exclusive
        if 'TSV_Groups' in self._tool_inputs:
            basename = super()._get_basename('TSV_Groups')
        elif 'TSV_Shared' in self._tool_inputs:
            basename = super()._get_basename('TSV_Shared')
        elif 'TSV_Counts' in self._tool_inputs:
            basename = super()._get_basename('TSV_Counts')
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(basename + '.count.summary')]
