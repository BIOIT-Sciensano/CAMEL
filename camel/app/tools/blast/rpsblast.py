from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.blast.blast import Blast


class RpsBlast(Blast):
    """
    Protein domain search
    """

    def __init__(self) -> None:
        """
        Initialize tool
            :return: None
        """
        super().__init__('rpsblast', '2.14.0')

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'DB_BLAST' not in self._tool_inputs:
            raise InvalidToolInputError('RPS BLAST Database (key: DB_BLAST) is required!')
        super()._check_input()
