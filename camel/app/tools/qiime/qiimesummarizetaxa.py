from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.qiime.qiime import Qiime


class QiimeSummarizeTaxa(Qiime):
    """
    The summarize_taxa.py script provides summary information of the representation of taxonomic groups within
    each sample.
    """

    def __init__(self):
        """
        Initialize tool
        :return: None
        """
        super().__init__('qiime_summarize_taxa', '1.9.1')

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - BIOM key is required
        - No additional keys allowed
        - Only two files allowed per key
        :return: None
        """
        if 'BIOM' not in self._tool_inputs:
            raise InvalidToolInputError('Invalid input files (keys) given for summarize_taxa: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['BIOM']) != 1:
            raise InvalidToolInputError('Invalid number (!= 1) of files in each key given for \
                                                 summarize_taxa: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = super(QiimeSummarizeTaxa, self)._get_basename('BIOM')
        self._tool_outputs['BIOM'] = []
        self._tool_outputs['TSV_Taxonomy'] = []
        for level in self.__get_levels():
            self._tool_outputs['BIOM'] += [ToolIOFile(basename + '_L' + str(level) + '.biom')]
            self._tool_outputs['TSV_Taxonomy'] += [ToolIOFile(basename + '_L' + str(level) + '.txt')]

    def __get_levels(self):
        """
        Returns the levels that are specified in the options or returns the default set
        :return: List of levels
        """
        if 'level' in self._parameters:
            return self._parameters['level'].value.split(',')
        else:
            return [2, 3, 4, 5, 6]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = '-i {}'.format(self._tool_inputs['BIOM'][0])
        input_string += ' -o {}'.format(self._folder)
        return input_string
