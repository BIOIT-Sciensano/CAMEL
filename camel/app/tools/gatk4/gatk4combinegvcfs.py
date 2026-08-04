from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4CombineGVCFs(GATK4):
    """
    =============================
    GATK CombineGVCFs 4.1.9.0
    =============================
    Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations

    Required inputs:
    ----------------
    'gVCF':             ToolIOFile object. Multiple gVCF files.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'gVCF'              ToolIOFile object. A single (merged) gVCF file
    """

    def __init__(self) -> None:
        """
        Initialize the CombineGVCFs tool
        :return: None
        """
        super().__init__('gatk4 CombineGVCFs', '4.1.9.0')

        self._required_inputs = ['gVCF', 'FASTA_REF']
        self._output_type = 'gVCF'

    def _set_input(self) -> None:
        """
        Set the input specification
        Overrides method in parent class.
        :return: None
        """
        super()._set_input()
        for f in self._tool_inputs['gVCF']:
            self._input_string += f"--variant {f.path} "
