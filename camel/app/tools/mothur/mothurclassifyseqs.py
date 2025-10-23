import os.path

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.mothur.mothur import Mothur


class MothurClassifySeqs(Mothur):
    """
    The classify.seqs command allows the user to use several different
    methods to assign their sequences to the taxonomy outline of their
    choice.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('mothur_classify_seqs', '1.39.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA, FASTA_Ref and TSV_Taxonomy keys are required
        - Additional allowed keys: 'TSV_Counts', 'TSV_Names', 'TSV_Groups'
        - Only one input file per key allowed
        - TSV_Counts, TSV_Names and TSV_Groups are mutually exclusive but this is not checked
        :return: None
        """
        super()._check_input()
        if not all(key in self._tool_inputs for key in ['FASTA', 'FASTA_Ref', 'TSV_Taxonomy']):
            raise InvalidToolInputError('Missing input files (keys) for Mothur '
                                                 'classify.seqs: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTA_Ref', 'TSV_Taxonomy', 'TSV_Counts', 'TSV_Names', 'TSV_Groups']:
                raise InvalidToolInputError('Invalid input key given for Mothur '
                                                     'classify.seqs: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidToolInputError('Invalid number (max = 1) of files given for Mothur \
                                                     classify.seqs: {!r}'.format(self._tool_inputs))

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = ['fasta={}'.format(self._tool_inputs['FASTA'][0]),
                 'reference={}'.format(self._tool_inputs['FASTA_Ref'][0]),
                 'taxonomy={}'.format(self._tool_inputs['TSV_Taxonomy'][0])]
        # TSV_Counts, TSV_Names and TSV_Groups are mutually exclusive
        if 'TSV_Counts' in self._tool_inputs:
            items.append('count={}'.format(self._tool_inputs['TSV_Counts'][0]))
        elif 'TSV_Names' in self._tool_inputs:
            items.append('name={}'.format(self._tool_inputs['TSV_Names'][0]))
        elif 'TSV_Groups' in self._tool_inputs:
            items.append('group={}'.format(self._tool_inputs['TSV_Groups'][0]))
        items.append('outputdir={}'.format(self._folder))
        return ', '.join(items)

    def _set_output(self):
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = super()._get_basename()
        # File name depends on the specified method option
        method_extension = self.__get_method_extension()
        tax_extension = self.__get_tax_extension()
        self._tool_outputs['TSV_Taxonomy'] = [ToolIOFile(basename + tax_extension + method_extension + '.taxonomy')]
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(basename + tax_extension + method_extension + '.tax.summary')]

    def __get_method_extension(self):
        """
        Checks whether a different method is specified in the options as
        the ouput file names are based on the method specified.
        Wang is the default method so '.wang' is returned if the knn option
        was not set.
        :return: String with extension
        """
        for name, parameter in self._parameters.items():
            if name == 'method' and parameter.value == 'knn':
                return '.knn'
        return '.wang'

    def __get_tax_extension(self):
        """
        Part of the taxonomy file name is used in the output. More specifically
        the part between the last two '.' is used or before the last '.' if
        there is no second '.' This method returns the relevant part of the
        taxonomy file name.
        :return: String with taxonomy extension used for output file naming
        """
        taxonomy = os.path.basename(self._tool_inputs['TSV_Taxonomy'][0].path)
        parts = taxonomy.split('.')
        return '.' + str(parts[-2])
