from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurListSeqs(Mothur):
    """
    The list.seqs command will write out the names of the sequences found within a fastq, fasta, name, group,
    count, list, or align.report file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_list_seqs')

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - One allowed key from: 'FASTA', 'FASTQ', 'TSV_Names', 'TSV_Counts', 'TSV_List', 'TSV_Groups'
        - Only one input file per key allowed
        :return: None
        """
        allowed_keys = {
            'FASTA',
            'FASTQ',
            'TSV_Names',
            'TSV_Counts',
            'TSV_List',
            'TSV_Groups',
        }
        counts = [
            True if key in allowed_keys else False for key in self._tool_inputs.keys()
        ].count(True)
        if counts > 1:
            raise InvalidToolInputError(
                f'Too many input keys given for Mothur list.seqs: {self._tool_inputs}'
            )
        if counts == 0:
            raise InvalidToolInputError(
                f'Invalid input key given for Mothur list.seqs: {self._tool_inputs}'
            )
        self._required_input = [
            key for key in self._tool_inputs.keys() if key in allowed_keys
        ]
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = []
        input_parameters = {
            'FASTA': 'fasta=',
            'FASTQ': 'fastq=',
            'TSV_Names': 'name=',
            'TSV_Counts': 'count=',
            'TSV_Groups': 'group=',
            'TSV_List': 'list=',
        }
        # Based on the key the correct option flag is added to the input string
        for key, input_files in self._tool_inputs.items():
            items.append(f'{input_parameters[key]}{input_files[0].path}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output file, and fills the common stream object with them
        :return: None
        """
        output_extension = '.accnos'
        input_key = next(iter(self._tool_inputs))
        self._tool_outputs['TSV_Accnos'] = [
            ToolIOFile(self._get_basename(input_key).with_suffix(output_extension))
        ]
