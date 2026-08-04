from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurRemoveRare(Mothur):
    """
    The remove.rare command removes OTUs at a specified rarity (number of observations in the dataset) and outputs a
    new file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_remove_rare')
        self._required_input = ['TSV_List']
        self._optional_input = ['TSV_Counts']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        items = [f"list={self._tool_inputs['TSV_List'][0]}", f"outputdir={self._folder}"]
        if 'TSV_Counts' in self._tool_inputs:
            items.append(f"count={self._tool_inputs['TSV_Counts'][0]}")
        return ', '.join(items)

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        labels = self._get_labels()
        basename = self._get_basename('TSV_List')
        self._tool_outputs['TSV_List'] = []
        # Only the first label in the file is used in case a list file is given as input
        self._tool_outputs['TSV_List'].append(ToolIOFile(basename.with_suffix(f'.{labels[0]}.pick.list')))
        if 'TSV_Counts' in self._tool_inputs:
            self._tool_outputs['TSV_Counts'] = [ToolIOFile(self._get_basename('TSV_Counts').with_suffix('.pick.count_table'))]
