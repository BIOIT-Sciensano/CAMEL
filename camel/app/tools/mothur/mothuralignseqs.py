from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurAlignSeqs(Mothur):
    """
    The align.seqs command aligns a user-supplied fasta-formatted candidate sequence file to a user-supplied
    fasta-formatted template alignment. The general approach is to i) find the closest template for each candidate
    using kmer searching, blastn, or suffix tree searching; ii) to make a pairwise alignment between the candidate and
    de-gapped template sequences using the Needleman-Wunsch, Gotoh, or blastn algorithms; and iii) to re-insert gaps
    to the candidate and template pairwise alignments using the NAST algorithm so that the candidate sequence alignment
    is compatible with the original template alignment.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_align_seqs')
        self._required_input = ['FASTA', 'FASTA_Ref']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        return ', '.join(
            [
                f"fasta={self._tool_inputs['FASTA'][0]}",
                f"reference={self._tool_inputs['FASTA_Ref'][0]}",
                f"outputdir={self._folder}",
            ]
        )

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        basename = self._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix('.align'))]
        self._tool_outputs['TSV_Report'] = [
            ToolIOFile(basename.with_suffix('.align_report'))
        ]
