from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurDistSeqs(Mothur):
    """
    The dist.seqs command will calculate uncorrected pairwise distances between aligned DNA sequences. This approach
    is better than the commonly used DNADIST because the distances are not stored in RAM, rather they are printed
    directly to a file. Furthermore, it is possible to ignore "large" distances that one might not be interested in.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_dist_seqs')
        self._required_input = ['FASTA']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}", f'outputdir={self._folder}']
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename()
        self._tool_outputs['DIST'] = [ToolIOFile(basename.with_suffix('.dist'))]
