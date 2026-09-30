from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.gatk.gatk import GATK


class GATKIndelRealigner(GATK):
    """
    Performs local realignment of reads around indels.
    Note: only required for variant callers without a haplotype assembly step (e.g. LoFreq, MuTect1), not for
    HaplotypeCaller or MuTect2.

    Required inputs:
    - BAM: Alignment file
    - FASTA_REF: Reference genome
    - TXT_realign_intervals: Target intervals created by RealignerTargetCreator

    Optional inputs:
    - VCF_KNOWN: Known indels
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - BAM: Realigned alignment file
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk IndelRealigner')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM', 'FASTA_REF', 'TXT_realign_intervals'])
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [
            *super()._get_input_options(),
            f"-I {self._tool_inputs['BAM'][0].path}",
            f"-targetIntervals {self._tool_inputs['TXT_realign_intervals'][0].path}"
        ]
        if 'VCF_KNOWN' in self._tool_inputs:
            options.append(f"-known {self._tool_inputs['VCF_KNOWN'][0].path}")
        return options

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self.get_param_value('output'))]
