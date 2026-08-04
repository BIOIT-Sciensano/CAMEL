from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurMakeShared(Mothur):
    """
    The make.shared command reads a list and group file or biom file and creates a .shared file as well as a rabund
    file for each group.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_make_shared')
        self._required_input = []
        self._optional_input = ['TSV_Groups', 'TSV_Counts']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - Either TSV_List or BIOM is required (not both)
        - Additional allowed keys are: 'TSV_Groups', 'TSV_Counts'
        - Only one input file per key allowed
        :return: None
        """
        if 'TSV_List' in self._tool_inputs:
            self._required_input = ['TSV_List']
        elif 'BIOM' in self._tool_inputs:
            self._required_input = ['BIOM']
        else:
            raise InvalidToolInputError('Invalid input files given. Please provide either TSV_List or BIOM.')
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'TSV_List' in self._tool_inputs:
            items.append(f"list={self._tool_inputs['TSV_List'][0]}")
        elif 'BIOM' in self._tool_inputs:
            items.append(f"biom={self._tool_inputs['BIOM'][0]}")
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them.
        REMARK: According to the documentation more output files will be created when
        a group file is given but it is not documented which ones.
        :return: None
        """
        # Either TSV_List or BIOM is given
        basename = self._get_basename('TSV_List') if 'TSV_List' in self._tool_inputs else self._get_basename('BIOM')
        self._tool_outputs['TSV_Shared'] = [ToolIOFile(basename.with_suffix('.shared'))]
