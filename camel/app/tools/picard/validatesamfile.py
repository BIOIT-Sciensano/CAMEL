from camel.app.tools.picard.picard import Picard


class ValidateSamFile(Picard):
    """
    Class for Picard ValidateSamFile function
    """

    def __init__(self):
        """
        Initialize a picard tool
                :return: None
        """
        super().__init__('Picard ValidateSamFile', '2.23.3')
        self._output_type = 'TXT_report'
