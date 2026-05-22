from camelcore.app.utils import fastqutils

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class FastqStats(Tool):
    """
    Calculates statistics for input FASTQ files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('FastqStats', '0.1')

    def _check_input(self) -> None:
        """
        Check if the provided tool input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError("FASTQ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._informs['stats'] = []
        for fastq_file in self._tool_inputs['FASTQ']:
            self._informs['stats'].append({
                'path': str(fastq_file.path),
                'nb_of_sequences': fastqutils.count_reads(fastq_file.path),
                'nb_of_bases': fastqutils.count_bases(fastq_file.path)
            })
