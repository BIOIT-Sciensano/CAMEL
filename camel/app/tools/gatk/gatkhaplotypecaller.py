from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.gatk.gatk import GATK


class GATKHaplotypeCaller(GATK):
    """
    Class for GATK HaplotypeCaller function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('gatk HaplotypeCaller', '3.7')

        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'VCF'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKHaplotypeCaller, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += f"-I {bam_file} "

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        super(GATKHaplotypeCaller, self)._set_output()

        if 'bamOutput' in self._parameters:
            self._tool_outputs['BAM'] = [ToolIOFile(self._folder / self.get_param_value('bamOutput'))]
