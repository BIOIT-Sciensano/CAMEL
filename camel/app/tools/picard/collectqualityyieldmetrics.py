from camel.app.tools.picard.picard import Picard


class CollectQualityYieldMetrics(Picard):
    """
    Class for picard CollectQualityYieldMetrics function to calculate a set of metrics used to describe the
    general quality of a BAM file
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard CollectQualityYieldMetrics', '2.23.3')
        self._output_type = 'TXT'
