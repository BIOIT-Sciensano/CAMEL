from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.gatk.gatk import GATK


class GATKBaseRecalibrator(GATK):
    """
    Detects systematic errors in base quality scores and generates a recalibration table.

    Required inputs:
    - BAM: Alignment file
    - FASTA_REF: Reference genome
    - VCF_KNOWN_SNPS and/or VCF_KNOWN_INDELS: Known variant sites

    Optional inputs:
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - TXT_RecalibrationTable: Recalibration table
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk BaseRecalibrator')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM', 'FASTA_REF'])
        if not any(key in self._tool_inputs for key in ('VCF_KNOWN_SNPS', 'VCF_KNOWN_INDELS')):
            raise InvalidToolInputError("Known sites input is required (VCF_KNOWN_SNPS and/or VCF_KNOWN_INDELS)")
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [*super()._get_input_options(), f"-I {self._tool_inputs['BAM'][0].path}"]
        for key in ('VCF_KNOWN_SNPS', 'VCF_KNOWN_INDELS'):
            options.extend(f'-knownSites {vcf_file.path}' for vcf_file in self._tool_inputs.get(key, []))
        return options

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['TXT_RecalibrationTable'] = [
            ToolIOFile(self.folder / self.get_param_value('recal_table_output'))]
