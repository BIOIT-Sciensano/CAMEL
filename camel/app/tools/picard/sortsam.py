from camel.app.tools.picard.picard import Picard


class SortSam(Picard):
    """
    Class for Picard SortSam function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard SortSam', '2.23.3')
