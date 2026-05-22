from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurCountSeqs(Mothur):
    """
    The count.seqs command counts the number of sequences represented
    by the representative sequence in a name file. If a group file is
    given, it will also provide the group count breakdown.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_count_seqs', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Names key is required
        - Only one file for TSV_Names is allowed
        - Only one addtional key is allowed: TSV_Groups
        :return: None
        """
        super()._check_input()
        if 'TSV_Names' not in self._tool_inputs:
            raise InvalidToolInputError('No input file given for Mothur count.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['TSV_Names']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                 count.seqs: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError('Too many input keys given for Mothur count.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['TSV_Names', 'TSV_Groups']:
                raise InvalidToolInputError('Invalid input key given for Mothur count.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/,
        outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['name={}'.format(self._tool_inputs['TSV_Names'][0])]
        # Only TSV_Groups can be an additional input key
        if 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename('TSV_Names')
        self._tool_outputs['TSV_Counts'] = [ToolIOFile(basename + '.count_table')]
