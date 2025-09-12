from camel.app.tools.blast.blast import Blast


class Blastp(Blast):
    """
    BlastP compares a protein query to a protein database.
    """

    def __init__(self) -> None:
        """
        Initialize tool.
        :return: None
        """
        super().__init__('blastp', '2.14.0')
