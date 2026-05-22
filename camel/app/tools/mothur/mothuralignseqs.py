from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
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

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_align_seqs', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA and FASTA_REF keys are required
        - Only one input file per key allowed
        - No other input keys are allowed
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs or 'FASTA_Ref' not in self._tool_inputs:
            raise InvalidToolInputError(f'Not enough valid input files given for {self.name}: {self._tool_inputs}')
        if len(self._tool_inputs['FASTA']) != 1 or len(self._tool_inputs['FASTA_Ref']) != 1:
            raise InvalidToolInputError(f'Invalid number (max = 1) of files per key given for {self.name}: {self._tool_inputs}')
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError(f'Too many input keys given for {self.name}: {self._tool_inputs}')

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f'fasta={self._tool_inputs["FASTA"][0]}',
                 f'reference={self._tool_inputs["FASTA_Ref"][0]}',
                 f'outputdir={self._folder}']
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the output file object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.align')]
        self._tool_outputs['TSV_Report'] = [ToolIOFile(basename + '.align.report')]
