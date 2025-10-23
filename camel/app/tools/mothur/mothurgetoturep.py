from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurGetOturep(Mothur):
    """
    The get.oturep command generates a fasta-formatted sequence file containing only a representative sequence
    for each OTU.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_get_oturep', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - PHY or DIST key is required
        - TSV_List key is required
        - In case of DIST key as input, TSV_Counts is also required
        - Additional allowed keys: 'TSV_Counts', 'FASTA', 'TSV_Names', 'TSV_Groups'
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        if not (('PHY' in self._tool_inputs) != ('DIST' in self._tool_inputs)) \
                or 'TSV_List' not in self._tool_inputs:
            raise InvalidToolInputError('Missing input files (keys) for Mothur '
                                                 'get.oturep: {!r}'.format(self._tool_inputs))
        if 'TSV_Counts' in self._tool_inputs and 'TSV_Names' in self._tool_inputs:
            raise InvalidToolInputError('Both Count and Names file is not allowed for Mothur '
                                                 'get.oturep: {!r}'.format(self._tool_inputs))
        if 'DIST' in self._tool_inputs and ('TSV_Counts' not in self._tool_inputs and 'TSV_Names' not in self._tool_inputs):
            raise InvalidToolInputError('TSV_Counts or TSV_Names required of column distance matrix is given '
                                                 'for Mothur get.oturep: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['PHY', 'DIST', 'FASTA', 'TSV_Counts', 'TSV_Names', 'TSV_Groups', 'TSV_List']:
                raise InvalidToolInputError('Invalid input key given for Mothur get.oturep: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     get.oturep: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        if 'PHY' in self._tool_inputs:
            items.append('phylip={}'.format(self._tool_inputs['PHY'][0]))
        else:
            items.append('column={}'.format(self._tool_inputs['DIST'][0]))
        items.append('list={}'.format(self._tool_inputs['TSV_List'][0]))
        if 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        elif 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        if 'FASTA' in self._tool_inputs:
            items.append('fasta={}'.format(self._tool_inputs['FASTA'][0]))
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename('TSV_List')
        labels = super()._get_labels()
        self._tool_outputs['FASTA'] = []
        for label in labels:
            self._tool_outputs['FASTA'] += [ToolIOFile(basename + '.' + label + '.rep.fasta')]
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = []
            for label in labels:
                self._tool_outputs['TSV_Counts'] += [ToolIOFile(basename + '.' + label + '.rep.count_table')]
        elif 'TSV_Names' in self._tool_inputs:
            self._tool_outputs['TSV_Names'] = []
            for label in labels:
                self._tool_outputs['TSV_Names'] += [ToolIOFile(basename + '.' + label + '.rep.names')]
