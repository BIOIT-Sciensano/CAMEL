from camel.app.core.utils import toolutils
from camel.app.core.tool import Tool


class SeqtkSize(Tool):
    """
    Reports the number sequences and bases in a FASTQ file.
    """

    def __init__(self) -> None:
        """
        Initialize seqtk subsample
        :return: None
        """
        super().__init__('Seqtk size', '1.4')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['FASTQ'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Function to run seqtk subsample
        :return: None
        """
        self._informs['stats'] = []
        for path_fq in [x.path for x in self._tool_inputs['FASTQ']]:
            self._command.command = ' '.join([self._tool_command, str(path_fq)])
            self._execute_command()
            nb_reads, nb_bases = self._command.stdout.strip().split('\t')
            self.informs['stats'].append({
                'path': str(path_fq),
                'nb_of_bases': int(nb_bases),
                'nb_of_sequences': int(nb_reads)
            })
