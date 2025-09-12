from camel.app.tools.gatk.gatk import GATK


class GATKGenotypeGVCFs(GATK):
    """
    Class for GATK GenotypeGVCFs function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('gatk GenotypeGVCFs', '3.7')

        self._required_inputs = ['gVCF', 'FASTA_REF']
        self._output_type = 'VCF_MultipleSample'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKGenotypeGVCFs, self)._set_input()

        for f in self._tool_inputs['gVCF']:
            self._input_string += "--variant {} ".format(f.path)
