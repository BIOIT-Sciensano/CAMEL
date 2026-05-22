import os.path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.qiime.qiime import Qiime


class QiimePickRepSet(Qiime):
    """
    After picking OTUs, you can then pick a representative set of sequences. For each OTU, you will end up with one
    sequence that can be used in subsequent analyses.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_pick_rep_set', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Otu key is required
        - FASTA is allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'TSV_Otu' not in self._tool_inputs:
            raise InvalidToolInputError(f'Invalid input files (keys) given for pick_rep_set: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Otu']:
                raise InvalidToolInputError(f'Invalid input key given for pick_rep_set: {self._tool_inputs!r}')
            if len(self._tool_inputs[key]) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files in each key given for \
                                                     pick_rep_set: {self._tool_inputs!r}')

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = super(QiimePickRepSet, self)._get_basename('TSV_Otu')
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / f'{basename}_rep_set.fasta')]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        basename = super(QiimePickRepSet, self)._get_basename('TSV_Otu')
        input_string = '-i {}'.format(self._tool_inputs['TSV_Otu'][0])
        if 'FASTA' in self._tool_inputs:
            input_string += ' -f {}'.format(self._tool_inputs['FASTA'][0])
        input_string += ' -o {}'.format(os.path.join(self._folder, basename + '_rep_set.fasta'))
        return input_string
