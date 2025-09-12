from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveLineage(Mothur):
    """
    The remove.lineage command reads a taxonomy file and a taxon and generates a new file that contains only the
    sequences not containing that taxon.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_remove_lineage', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV_Taxonomy key is required
        - Other allowed keys are: 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
          'TSV_AlignReport', 'TSV_List'
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        if 'TSV_Taxonomy' not in self._tool_inputs:
            raise InvalidToolInputError('Not enough valid input files given for Mothur '
                                                 'remove.lineage: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
                           'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy']:
                raise InvalidToolInputError('Invalid input key given for Mothur remove.lineage: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     remove.lineage: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['taxonomy={}'.format(self._tool_inputs['TSV_Taxonomy'][0])]
        input_parameters = {'FASTA': 'fasta=',
                            'TSV_Names': 'name=',
                            'TSV_Counts': 'count=',
                            'TSV_Groups': 'group=',
                            'TSV_AlignReport': 'alignreport=',
                            'TSV_List': 'list='}
        for key, input_files in self._tool_inputs.items():
            # Based on the key the correct option flag is added to the input string
            if key != 'TSV_Taxonomy':
                items.append('{}{}'.format(input_parameters[key], input_files[0].path))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # The extension of the output files is not always built in the same way. A dictionary
        # holds the suffixes that can be used when getting the base name
        output_extensions = {'FASTA': ['.', '.pick.fasta'],
                             'TSV_Names': ['.', '.pick.names'],
                             'TSV_Counts': ['.', '.pick.count_table'],
                             'TSV_Groups': ['.', '.pick.groups'],
                             'TSV_AlignReport': ['.align.report', '.pick.align.report'],
                             'TSV_List': ['.', '.pick.list'],
                             'TSV_Taxonomy': ['.', '.pick.taxonomy']}
        for key, input_files in self._tool_inputs.items():
            basename = super()._get_basename(key, output_extensions[key][0])
            self._tool_outputs[key] = [ToolIOFile(basename + output_extensions[key][1])]
