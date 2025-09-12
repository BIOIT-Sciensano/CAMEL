from camel.app.tools.picard.picard import Picard


class CollectWgsMetrics(Picard):
    """
    Class for picard CollectWgsMetrics function. Collect metrics about coverage and performance of WGS experiments.
    """

    def __init__(self):
        """
        Initialize a picard tool
                :return: None
        """
        super().__init__('Picard CollectWgsMetrics', '2.23.3')
        self._required_inputs = ['BAM', 'SAM', 'FASTA_REF']
        self._output_type = 'TXT_metrics'

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        super()._set_input()
        if 'COVERAGE_INTERVALS' in self._tool_inputs:
            self._input_string += f" INTERVALS={self._tool_inputs['COVERAGE_INTERVALS'][0].path}"
