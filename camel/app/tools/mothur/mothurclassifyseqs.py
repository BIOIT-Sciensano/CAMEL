from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurClassifySeqs(Mothur):
    """
    The classify.seqs command allows the user to use several different
    methods to assign their sequences to the taxonomy outline of their
    choice.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_classify_seqs')
        self._required_input = ['FASTA', 'FASTA_Ref', 'TSV_Taxonomy']
        self._optional_input = ['TSV_Counts', 'TSV_Names', 'TSV_Groups']

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid cfr. super
        Additionally:
        - TSV_Counts, TSV_Names and TSV_Groups are mutually exclusive
        :return: None
        """
        super()._check_input()
        if len([x for x in self._tool_inputs if x in self._optional_input]) != 1:
            raise InvalidToolInputError(
                "Invalid input files (keys) given. 'TSV_Counts', 'TSV_Names' and 'TSV_Groups' are mutually exclusive.")

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}",
                 f"reference={self._tool_inputs['FASTA_Ref'][0]}",
                 f"taxonomy={self._tool_inputs['TSV_Taxonomy'][0]}"]
        # TSV_Counts, TSV_Names and TSV_Groups are mutually exclusive
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        elif 'TSV_Names' in self._tool_inputs:
            items.append(f"name={self._tool_inputs['TSV_Names'][0]}")
        elif 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename()
        # File name depends on the specified method option
        method_extension = self.__get_method_extension()
        tax_extension = self.__get_tax_extension()
        self._tool_outputs['TSV_Taxonomy'] = [ToolIOFile(basename.with_suffix(f'{tax_extension}{method_extension}.taxonomy'))]
        self._tool_outputs['TSV_Summary'] = [ToolIOFile(basename.with_suffix(f'{tax_extension}{method_extension}.tax.summary'))]

    def __get_method_extension(self) -> str:
        """
        Checks whether a different method is specified in the options as
        the ouput file names are based on the method specified.
        Wang is the default method so '.wang' is returned if the knn option
        was not set.
        :return: String with extension
        """
        return '.knn' if self._parameters.get('method', 'wang') == 'knn' else '.wang'

    def __get_tax_extension(self) -> str:
        """
        Part of the taxonomy file name is used in the output. More specifically
        the part between the last two '.' is used or before the last '.' if
        there is no second '.' This method returns the relevant part of the
        taxonomy file name.
        :return: String with taxonomy extension used for output file naming
        """
        taxonomy = str(self._tool_inputs['TSV_Taxonomy'][0].basename)
        parts = taxonomy.split('.')
        return '.' + str(parts[-2])
