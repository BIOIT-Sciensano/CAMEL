from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurDistSeqs(Mothur):
    """
    The dist.seqs command will calculate uncorrected pairwise distances between aligned DNA sequences. This approach
    is better than the commonly used DNADIST because the distances are not stored in RAM, rather they are printed
    directly to a file. Furthermore, it is possible to ignore "large" distances that one might not be interested in.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_dist_seqs', version=None)

    def check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - No additional keys are allowed
        - Only one FASTA file is allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for Mothur '
                                                 'dist.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Invalid input keys given for Mothur dist.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files in each key given for Mothur '
                                                 'dist.seqs: {!r}'.format(self._tool_inputs))

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
        self._tool_outputs['DIST'] = [ToolIOFile(basename + '.dist')]
