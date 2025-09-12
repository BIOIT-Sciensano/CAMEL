import pandas as pd

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqkit.seqkitbase import SeqkitBase


class SeqkitStats(SeqkitBase):
    """
    Seqkit stats reports simple statistics of FASTA/Q files.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self) -> None:
        """
        Initializes this tool.
                :return: None
        """
        super().__init__('Seqkit stats', version=None)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in SeqkitStats.INPUT_KEYS):
            raise InvalidToolInputError(f'{" or ".join(SeqkitStats.INPUT_KEYS)} input is required.')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        tsv_out = self.folder / self._parameters['output_filename'].value
        self.__build_command()
        self._execute_command()
        data = pd.read_table(tsv_out)
        for k, v in next(iter(data.to_dict('records'))).items():
            self._informs[str(k)] = v
        self._tool_outputs['TSV'] = [ToolIOFile(tsv_out)]

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        self._command.command = ' '.join([
            self._tool_command,
            '--all', '--tabular',
            *([str(f.path) for f in self._tool_inputs[input_key]]),
            *self._build_options()
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
