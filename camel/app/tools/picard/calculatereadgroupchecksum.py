from camel.app.tools.picard.picard import Picard


class CalculateReadGroupChecksum(Picard):
    """
    Class for Picard CalculateReadGroupChecksum 2.23.3
    Creates a hash code based on the read groups (RG)
    Output: 'TXT_checksum'  ToolIOFile object. File to which the hash code will be written
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard CalculateReadGroupChecksum', '2.23.3')
        self._output_type = 'TXT_checksum'
