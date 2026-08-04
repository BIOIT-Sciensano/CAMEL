from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurPcrSeqs(Mothur):
    """
    The pcr.seqs will trim inputted sequences based on a variety of user-defined options.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_pcr_seqs')
        self._required_input = ['FASTA']
        self._optional_input = ['TSV_Oligos']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        items = [f"fasta={self._tool_inputs['FASTA'][0]}"]
        if 'TSV_Oligos' in self._tool_inputs:
            items.append(f"oligos={self._tool_inputs['TSV_Oligos'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        basename = self._get_basename()
        extension = self._tool_inputs['FASTA'][0].file_extension
        self._tool_outputs['FASTA'] = [ToolIOFile(basename.with_suffix(f'.pcr{extension}'))]
        if basename.with_suffix('.bad.accnos').exists():
            self._tool_outputs['TEXT'] = [ToolIOFile(basename.with_suffix('.bad.accnos'))]
