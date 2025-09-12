from camel.app.tools.picard.picard import Picard


class MergeVCFs(Picard):
    """
    Class for picard MergeVCFs function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard MergeVCFs', '2.23.3')

        self._required_inputs = ['VCF']
        self._output_type = 'VCF'

    def _set_input(self) -> None:
        """
        Set the input specification. This method handles on or more VCF files
        :return: None
        """
        self._input_string = "".join(f'I={vcf.path} ' for vcf in self._tool_inputs["VCF"])
