import os.path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.qiime.qiime import Qiime


class QiimeMakeOtuTable(Qiime):
    """
    The script make_otu_table.py tabulates the number of times an OTU is found in each sample, and adds the taxonomic
    predictions for each OTU in the last column if a taxonomy file is supplied.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_make_otu_table', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Otu key is required
        - TSV_Taxonomy is allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'TSV_Otu' not in self._tool_inputs:
            raise InvalidToolInputError(f'Invalid input files (keys) given for make_otu_table: {self._tool_inputs!r}')
        for key, input_files in self._tool_inputs.items():
            if key not in ['TSV_Taxonomy', 'TSV_Otu']:
                raise InvalidToolInputError(f'Invalid input key given for make_otu_table: {self._tool_inputs!r}')
            if len(self._tool_inputs[key]) != 1:
                raise InvalidToolInputError(f'Invalid number (max = 1) of files in each key given for \
                                                     make_otu_table: {self._tool_inputs!r}')

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = super(QiimeMakeOtuTable, self)._get_basename('TSV_Otu')
        self._tool_outputs['BIOM'] = [ToolIOFile(self._folder / basename + '_otu_table.biom')]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        basename = super(QiimeMakeOtuTable, self)._get_basename('TSV_Otu')
        input_string = '-i {}'.format(self._tool_inputs['TSV_Otu'][0])
        if 'TSV_Taxonomy' in self._tool_inputs:
            input_string += ' -t {}'.format(self._tool_inputs['FASTA'][0])
        input_string += ' -o {}'.format(os.path.join(self._folder, basename + '_otu_table.biom'))
        return input_string
