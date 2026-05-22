from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurFilterSeqs(Mothur):
    """
    filter.seqs removes columns from alignments based on a criteria defined by the user.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_filter_seqs')
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
        self._tool_outputs['FASTA'] = [ToolIOFile(self._get_basename().with_suffix('.filter.fasta'))]
