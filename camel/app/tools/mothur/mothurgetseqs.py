from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurGetSeqs(Mothur):
    """
    The get.seqs command takes a list of sequence names (.accnos file) and either a fastq, fasta, name, group,
    list, count or align.report file to generate a new file that contains only the sequences in the list.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_get_seqs')
        self._required_input = ['TSV_Accnos']
        self._optional_input = [
            'FASTA',
            'TSV_Names',
            'TSV_Counts',
            'TSV_List',
            'TSV_Taxonomy',
            'TSV_Groups',
        ]

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - TSV_Accnos is required
        - One additional allowed key from: 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_List', 'TSV_Taxonomy', 'TSV_Groups'
        - Only one input file per key allowed
        :return: None
        """
        super()._check_input()
        allowed_keys = self._required_input + self._optional_input
        if [
            True if key in allowed_keys else False for key in self._tool_inputs.keys()
        ].count(True) > 2:
            raise InvalidToolInputError('Too many input keys given.')

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        input_parameters = {
            'TSV_Accnos': 'accnos=',
            'FASTA': 'fasta=',
            'TSV_Names': 'name=',
            'TSV_Counts': 'count=',
            'TSV_Groups': 'group=',
            'TSV_List': 'list=',
            'TSV_Taxonomy': 'taxonomy=',
        }
        # Based on the key the correct option flag is added to the input string
        for key, input_files in self._tool_inputs.items():
            items.append(f'{input_parameters[key]}{input_files[0].path}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        output_extensions = {
            'FASTA': '.pick.fasta',
            'TSV_Taxonomy': '.pick.taxonomy',
            'TSV_Names': '.pick.names',
            'TSV_Counts': '.pick.count_table',
            'TSV_Groups': '.pick.groups',
            'TSV_List': '.pick.list',
        }
        # Based on the key the correct output file is added to the input string
        for key, input_files in self._tool_inputs.items():
            if key == 'TSV_Accnos':
                continue
            basename = self._get_basename(key)
            self._tool_outputs[key] = [
                ToolIOFile(basename.with_suffix(output_extensions[key]))
            ]
