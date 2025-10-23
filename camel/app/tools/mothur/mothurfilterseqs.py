from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurFilterSeqs(Mothur):
    """
    filter.seqs removes columns from alignments based on a criteria defined by the user.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_filter_seqs', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only one FASTA file allowed
        - No other input keys allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('No valid input files given for Mothur filter.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                 filter.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Too many input keys given for Mothur filter.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0]),
                 'outputdir={}'.format(self._folder)]
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.filter.fasta')]
