from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurUniqueSeqs(Mothur):
    """
    The unique.seqs command returns only the unique sequences found in a
    fasta-formatted sequence file and a file that indicates those
    sequences that are identical to the reference sequence.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_unique_seqs')
        self._required_input = ['FASTA']
        self._optional_input = ['TSV_Counts']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix(f'.unique{self._get_extension()}'))]
        self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename.with_suffix('.count_table'))]
