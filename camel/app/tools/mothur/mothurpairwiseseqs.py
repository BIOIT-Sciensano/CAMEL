from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurPairwiseSeqs(Mothur):
    """
    The pairwise.seqs command will calculate uncorrected pairwise distances between sequences. The command will
    generate a column-formatted distance matrix that is compatible with the cluster command. The command is also able
    to generate a phylip-formatted distance matrix. There are several options for how to handle gap comparisons and
    terminal gaps.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_pairwise_seqs', '1.39.1')

    def _check_input(self):
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
                                                 'pairwise.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidToolInputError('Invalid input keys given for Mothur pairwise.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files in each key given for Mothur '
                                                 'pairwise.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={0}'.format(self._tool_inputs['FASTA'][0]),
                 'outputdir={0}'.format(self._folder)]
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['DIST'] = [ToolIOFile(basename + '.dist')]
