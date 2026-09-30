from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.gatk.gatk import GATK


class GATKPrintReads(GATK):
    """
    Writes the reads of an alignment file, e.g. to apply a base quality score recalibration table.

    Required inputs:
    - BAM: Alignment file
    - FASTA_REF: Reference genome

    Optional inputs:
    - BQSR: Recalibration table created by BaseRecalibrator
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - BAM: Alignment file
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk PrintReads')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM', 'FASTA_REF'])
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [*super()._get_input_options(), f"-I {self._tool_inputs['BAM'][0].path}"]
        if 'BQSR' in self._tool_inputs:
            options.append(f"-BQSR {self._tool_inputs['BQSR'][0].path}")
        return options

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self.get_param_value('bam_output'))]
