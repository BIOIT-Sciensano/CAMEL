from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.gatk.gatk import GATK


class GATKHaplotypeCaller(GATK):
    """
    Calls germline SNPs and indels via local re-assembly of haplotypes.

    Required inputs:
    - BAM: Alignment file
    - FASTA_REF: Reference genome

    Optional inputs:
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - VCF: Variants (or gVCF, depending on 'emit_ref_confidence')
    - BAM: Assembled haplotypes (only if 'bam_output' is set)
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk HaplotypeCaller')

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
        return [*super()._get_input_options(), f"-I {self._tool_inputs['BAM'][0].path}"]

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(self.folder / self.get_param_value('output'))]
        if self.get_param_value('bam_output') is not None:
            self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self.get_param_value('bam_output'))]
