from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurGetOturep(Mothur):
    """
    The get.oturep command generates a fasta-formatted sequence file containing only a representative sequence
    for each OTU.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_get_oturep')
        self._required_input = ['TSV_List']
        self._optional_input = ['TSV_Counts', 'FASTA', 'TSV_Names', 'TSV_Groups']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - PHY or DIST key is required (except if method=abundance)
        - TSV_List key is required
        - In case of DIST key as input, TSV_Counts or TSV_Names is also required
        - Additional allowed keys: 'TSV_Counts', 'FASTA', 'TSV_Names', 'TSV_Groups'
        - Only one input file per key allowed
        :return: None
        """
        if 'PHY' in self._tool_inputs:
            self._required_input.append('PHY')
        if 'DIST' in self._tool_inputs:
            self._required_input.append('DIST')
            if 'TSV_Counts' in self._tool_inputs:
                self._required_input.append('TSV_Counts')
            elif 'TSV_Names' in self._tool_inputs:
                self._required_input.append('TSV_Names')
            else:
                raise InvalidToolInputError('When using DIST, TSV_Counts or TSV_Names are also required.')
        # If method is abundance, PHY or DIST not required
        if not (('PHY' in self._tool_inputs) != ('DIST' in self._tool_inputs)):
            if not ('method' in self._parameters and self.get_param_value('method') == 'abundance'):
                raise InvalidToolInputError('Wrong input file combination (keys) for Mothur.')
        if 'TSV_Names' in self._tool_inputs and 'TSV_Counts' in self._tool_inputs:
            raise InvalidToolInputError("Invalid input files given. 'TSV_Names' and 'TSV_Counts' are mutually exclusive")
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'PHY' in self._tool_inputs:
            items.append(f"phylip={self._tool_inputs['PHY'][0]}")
        elif 'DIST' in self._tool_inputs:
            items.append(f"column={self._tool_inputs['DIST'][0]}")
        items.append(f"list={self._tool_inputs['TSV_List'][0]}")
        if 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        elif 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        if 'FASTA' in self._tool_inputs:
            items.append(f"fasta={self._tool_inputs['FASTA'][0]}")
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename('TSV_List')
        labels = self._get_labels()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix(f'.{label}.rep.fasta')) for label in labels]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename.with_suffix(f'.{label}.rep.count_table')) for label in labels]
        elif 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = [ToolIOFile(basename.with_suffix(f'.{label}.rep.names')) for label in labels]
