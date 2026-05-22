from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurRarefactionSingle(Mothur):
    """
    The rarefaction.single command will generate intra-sample rarefaction curves using a re-sampling without
    replacement approach.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_rarefaction_single', version=None)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_List key is required
        - Only one file for TSV_List is allowed
        - Use of .shared file not yet implemented
        :return: None
        """
        super()._check_input()
        if 'TSV_List' not in self._tool_inputs:
            raise InvalidToolInputError('No valid input file given for Mothur rarefaction.single: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['TSV_List']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                 rarefaction.single: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 1:
            raise InvalidToolInputError('Too many input keys given for Mothur rarefaction.single: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['list={}'.format(self._tool_inputs['TSV_List'][0]), 'outputdir={}'.format(self._folder)]
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename('TSV_List')
        self._tool_outputs['TSV'] = [ToolIOFile('{}.rarefaction'.format(basename))]
