from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurChimeraVsearch(Mothur):
    """
    The chimera.vsearch command reads a fasta file and reference file or a fasta and name or count file and outputs
    potentially chimeric sequences. The vsearch program is donated to the public domain,
    https://github.com/torognes/vsearch
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_chimera_vsearch', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Either TSV_Names or TSV_Counts is required
        - FASTA_Ref and TSV_Groups are allowed as additional input
        - Only one input file per key is allowed
        - The use of TSV_Names is not yet implemented (lack of documentation)
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur '
                                                 f'chimera.vsearch: {self._tool_inputs!r}')
        if 'TSV_Names' not in self._tool_inputs and 'TSV_Counts' not in self._tool_inputs:
            raise InvalidToolInputError('Missing input files (key) for Mothur '
                                                 f'chimera.vsearch: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Counts', 'TSV_Names', 'TSV_Groups', 'FASTA_Ref']:
                raise InvalidToolInputError('Invalid input key given for Mothur '
                                                     f'chimera.vsearch: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files in each key given for Mothur \
                                                     chimera.vsearch: {self._tool_inputs!r}')

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
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        if 'FASTA_Ref' in self._tool_inputs:
            items.append('reference={}'.format(self._tool_inputs['FASTA_Ref'][0]))
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['TSV_Chimeras'] = [ToolIOFile(basename + '.denovo.vsearch.chimeras')]
        self._tool_outputs['TSV_Accnos'] = [ToolIOFile(basename + '.denovo.vsearch.accnos')]
        if 'TSV_Names' in self._tool_inputs:
            raise RuntimeError('The use of a names file is not yet implemented for chimera.vsearch as the '
                               'outputs are unknown!')
