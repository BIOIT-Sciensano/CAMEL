from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveSeqs(Mothur):
    """
    The remove.seqs command takes a list of sequence names and either a fastq, fasta, name, group, list, count or
    align.report file to generate a new file that does not contain the sequences in the list.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_remove_seqs')
        self._required_input = ['TSV_Accnos']
        self._optional_input = [
            'FASTA',
            'TSV_Names',
            'TSV_Counts',
            'TSV_Groups',
            'TSV_AlignReport',
            'TSV_List',
            'TSV_Taxonomy',
            'TSV_Qfile',
            'FASTQ',
        ]

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - TSV_Accnos is required
        - In addition only two extra keys are permitted
        - Possible additional key comes from: 'FASTA', 'TSV_Names', 'TSV_Counts', 'TSV_Groups',
          'TSV_AlignReport', 'TSV_List', 'TSV_Taxonomy', 'TSV_Qfile', 'FASTQ'
        - Only one file per key is allowed
        :return: None
        """
        if len(self._tool_inputs.keys()) > 3:
            raise InvalidToolInputError(
                f'Too many input keys given for Mothur remove.seqs: {self._tool_inputs}'
            )
        super()._check_input()

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        # TSV_Accnos is required so will always be available
        items = [f"accnos={self._tool_inputs['TSV_Accnos'][0]}"]
        input_parameters = {
            'FASTA': 'fasta=',
            'FASTQ': 'fastq=',
            'TSV_Names': 'name=',
            'TSV_Counts': 'count=',
            'TSV_Groups': 'group=',
            'TSV_AlignReport': 'alignreport=',
            'TSV_List': 'list=',
            'TSV_Taxonomy': 'taxonomy=',
            'TSV_Qfile': 'qfile=',
        }
        for key, input_files in self._tool_inputs.items():
            # Only two keys are possible so the one that is not TSV_Accnos
            # will define the option flag that is needed
            if key != 'TSV_Accnos':
                items.append(f'{input_parameters[key]}{input_files[0].path}')
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        # The extension of the output files is not always built in the same way. A dictionary
        # holds the suffixes that can be used when getting the base name
        output_extensions = {
            'FASTA': [None, '.pick.fasta'],
            'FASTQ': [None, '.pick.fastq'],
            'TSV_Names': [None, '.pick.names'],
            'TSV_Counts': [None, '.pick.count_table'],
            'TSV_Groups': [None, '.pick.groups'],
            'TSV_AlignReport': [{'.report'}, '.pick.align.report'],
            'TSV_List': [None, '.pick.list'],
            'TSV_Taxonomy': [None, '.pick.taxonomy'],
            'TSV_Qfile': [None, '.pick.qual'],
        }
        for key, input_files in self._tool_inputs.items():
            # The TSV_Accnos file does not directly lead to output
            if key != 'TSV_Accnos':
                basename = self._get_basename(key, output_extensions[key][0])
                self._tool_outputs[key] = [
                    ToolIOFile(basename.with_suffix(output_extensions[key][1]))
                ]
