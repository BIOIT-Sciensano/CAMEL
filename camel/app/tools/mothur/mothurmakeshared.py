from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurMakeShared(Mothur):
    """
    The make.shared command reads a list and group file or biom file and creates a .shared file as well as a rabund
    file for each group.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_make_shared', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Either TSV_List or BIOM is required (not both)
        - Additional allowed keys are: 'TSV_Groups', 'TSV_Counts'
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        if not (('TSV_List' in self._tool_inputs) != ('BIOM' in self._tool_inputs)):
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur make.shared, only '
                                                 'TSV_List or BIOM allowed: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['TSV_List', 'BIOM', 'TSV_Groups', 'TSV_Counts']:
                raise InvalidToolInputError('Invalid input key given for Mothur make.shared: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     make.shared: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'TSV_List' in self._tool_inputs:
            items.append('list={}'.format(self._tool_inputs['TSV_List'][0]))
        elif 'BIOM' in self._tool_inputs:
            items.append('biom={}'.format(self._tool_inputs['BIOM'][0]))
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them.
        REMARK: According to the documentation more output files will be created when
        a group file is given but it is not documented which ones.
        :return: None
        """
        if 'TSV_List' in self._tool_inputs:
            basename = super()._get_basename('TSV_List')
        # Either TSV_List or BIOM is given
        else:
            basename = super()._get_basename('BIOM')
        self._tool_outputs['TSV_Shared'] = [ToolIOFile(basename + '.shared')]
