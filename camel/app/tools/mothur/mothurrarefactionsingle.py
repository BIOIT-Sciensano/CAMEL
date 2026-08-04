from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.mothur.mothur import Mothur


class MothurRarefactionSingle(Mothur):
    """
    The rarefaction.single command will generate intra-sample rarefaction curves using a re-sampling without
    replacement approach.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mothur_rarefaction_single')
        self._required_input = ['TSV_List']

    def _build_input_string(self) -> str:
        """
        Creates the string with the input files and input/output directories
        Example: fasta=File1.trim.contig.fasta, outputdir=/test/data/outputdir
        :return: String with the input parameters
        """
        return f"list={self._tool_inputs['TSV_List'][0]}, outputdir={self._folder}"

    def _set_output(self) -> None:
        """
        Sets the name of the output files, and fills the common stream object with them
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(self._get_basename('TSV_List').with_suffix('.rarefaction'))]
