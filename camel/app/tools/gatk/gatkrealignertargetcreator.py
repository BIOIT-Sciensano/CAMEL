from camel.app.loggers import logger
from camel.app.tools.gatk.gatk import GATK


class GATKRealignerTargetCreator(GATK):
    """
    Class for GATK RealignerTargetCreator function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('gatk RealignerTargetCreator', '3.7')

        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'TXT_realign_intervals'
        logger.info(""" NOTE: From GATK 3.6 on, Indel realignment is no longer necessary for variant discovery if you plan to use a variant
            caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still
            required when using legacy callers such as UnifiedGenotyper or the original MuTect. """)

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKRealignerTargetCreator, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += f"-I {bam_file} "

        # to work with only known regions from VCF_KNOWN
        if 'VCF_KNOWN' in self._tool_inputs:
            self._input_string += "-known {} ".format(self._tool_inputs['VCF_KNOWN'][0].path)
