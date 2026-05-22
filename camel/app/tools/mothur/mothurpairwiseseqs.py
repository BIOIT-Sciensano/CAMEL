from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurPairwiseSeqs(Mothur):
    """
    The pairwise.seqs command will calculate uncorrected pairwise distances between sequences. The command will
    generate a column-formatted distance matrix that is compatible with the cluster command. The command is also able
    to generate a phylip-formatted distance matrix. There are several options for how to handle gap comparisons and
    terminal gaps.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_pairwise_seqs')
        self._required_input = ['FASTA']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        return f"fasta={self._tool_inputs['FASTA'][0]}, outputdir={self._folder}"

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['DIST'] = [ToolIOFile(self._get_basename().with_suffix('.dist'))]
