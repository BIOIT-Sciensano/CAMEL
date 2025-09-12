from pathlib import Path

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveRare(Mothur):
    """
    The remove.rare command removes OTUs at a specified rarity (number of observations in the dataset) and outputs a
    new file.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_remove_rare', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_List key is required
        - Only one file for TSV_List is allowed
        - TSV_Counts is the only extra input key that is allowed
        - Rabund and sabund not yet implemented
        :return: None
        """
        super()._check_input()
        if 'TSV_List' not in self._tool_inputs:
            raise InvalidToolInputError('No valid input file given for Mothur remove.rare: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['TSV_List']) != 1:
            raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                 remove.rare: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidToolInputError('Too many input keys given for Mothur remove.rare: {!r}'.format(self._tool_inputs))
        for key in self._tool_inputs.keys():
            if key not in {'TSV_List', 'TSV_Counts'}:
                raise InvalidToolInputError('Invalid input key given for Mothur remove.rare: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = ['list={}'.format(self._tool_inputs['TSV_List'][0]), 'outputdir={}'.format(self._folder)]
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        labels = self._get_labels()
        basename = self._get_basename('TSV_List')
        self._tool_outputs['TSV_List'] = []
        # Only the first label in the file is used in case a list file is given as input
        self._tool_outputs['TSV_List'].append(ToolIOFile(Path(f'{basename}.{labels[0]}.pick.list')))
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(Path(f'{self._get_basename("TSV_Counts")}.pick.count_table'))]
