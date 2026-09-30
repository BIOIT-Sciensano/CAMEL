from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.gatk.gatk import GATK


class GATKGenotypeGVCFs(GATK):
    """
    Performs joint genotyping on one or more gVCF files produced by HaplotypeCaller.

    Required inputs:
    - gVCF: One or more gVCF files
    - FASTA_REF: Reference genome

    Optional inputs:
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - VCF_MultipleSample: Genotyped variants
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk GenotypeGVCFs')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['gVCF', 'FASTA_REF'])
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        return [
            *super()._get_input_options(),
            *[f'--variant {gvcf_file.path}' for gvcf_file in self._tool_inputs['gVCF']]
        ]

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['VCF_MultipleSample'] = [ToolIOFile(self.folder / self.get_param_value('output'))]
