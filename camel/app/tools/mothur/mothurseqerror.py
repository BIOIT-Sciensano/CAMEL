from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurSeqError(Mothur):
    """
    Seq.error measures the error rate
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_seq_error', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA and FASTA_REF are required
        - For now 'TSV_Names', 'TSV_Counts', and 'TSV_Qfile' are not implemented (lack of documentation)
        - Only one file per key allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs or 'FASTA_REF' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur '
                                                 f'seq.error: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key in ['TSV_Names', 'TSV_Counts', 'TSV_Qfile']:
                raise InvalidToolInputError('These input keys are not yet implemented as documentation for these keys'
                                                     f' is not yet available: {self._tool_inputs!r}')
            elif key not in ['FASTA', 'FASTA_REF']:
                raise InvalidToolInputError(f'Invalid input key given for Mothur seq.error: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files in each key given for Mothur \
                                                     seq.error: {self._tool_inputs!r}')

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0]),
                 'reference={}'.format(self._tool_inputs['FASTA_REF'][0]),
                 f'outputdir={self._folder}']
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs.update({
            'TSV_Summary': [ToolIOFile(basename + '.error.summary')],
            'TSV_Seq': [ToolIOFile(basename + '.error.seq')],
            'TSV_Chimera': [ToolIOFile(basename + '.error.chimera')],
            'TSV_SeqForward': [ToolIOFile(basename + '.error.seq.forward')],
            'TSV_SeqReverse': [ToolIOFile(basename + '.error.seq.reverse')],
            'TSV_Count': [ToolIOFile(basename + '.error.count')],
            'TSV_Matrix': [ToolIOFile(basename + '.error.matrix')]
        })
