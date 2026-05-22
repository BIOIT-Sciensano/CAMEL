from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurUniqueSeqs(Mothur):
    """
    The unique.seqs command returns only the unique sequences found in a
    fasta-formatted sequence file and a file that indicates those
    sequences that are identical to the reference sequence.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_unique_seqs', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Allowed keys are 'FASTA', 'TSV_Counts'
        - Only one input file allowed per key
        :return: None
        """
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError(f'No input file given for Mothur unique.seqs: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Counts']:
                raise InvalidToolInputError(f'Invalid input key given for Mothur unique.seqs: {self._tool_inputs!r}')
            if len(input_files) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files given for Mothur \
                                                     unique.seqs: {self._tool_inputs!r}')

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0])]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f'count={self._tool_inputs["TSV_Counts"][0]}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        self._tool_outputs['FASTA'] = [ToolIOFile(basename + '.unique' + self._get_extension())]
        # If no count table is given, a names file is created. Otherwise the count table is updated
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename + '.count_table')]
        else:
            self._tool_outputs['TSV_Names'] = [ToolIOFile(basename + '.names')]
