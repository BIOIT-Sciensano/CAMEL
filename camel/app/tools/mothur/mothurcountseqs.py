from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurCountSeqs(Mothur):
    """
    The count.seqs command counts the number of sequences represented
    by the representative sequence in a name file. If a group file is
    given, it will also provide the group count breakdown.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_count_seqs')
        self._required_input = ['TSV_Names']
        self._optional_input = ['TSV_Groups']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, inputdir=/test/data/input/, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = [f"name={self._tool_inputs['TSV_Names'][0]}"]
        # Only TSV_Groups can be an additional input key
        if 'TSV_Groups' in self._tool_inputs:
            items.append(f"group={self._tool_inputs['TSV_Groups'][0]}")
        items.append(f'outputdir={self._folder}')
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['TSV_Counts'] = [ToolIOFile(self._get_basename('TSV_Names').with_suffix('.count_table'))]
