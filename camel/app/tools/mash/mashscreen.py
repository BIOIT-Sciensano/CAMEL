from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class MashScreen(Tool):
    """
    Mash screen determines whether query sequences are within a larger mixture of sequences.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mash screen', '2.3')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['DB'], keys_allowed=['DB', 'FASTQ', 'FASTA'])
        if 'FASTQ' not in self._tool_inputs and 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTQ or FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        key_input = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        path_tsv_out = self.folder / 'mash_screen.tsv'
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['DB'][0].path),
            *(str(x.path) for x in self._tool_inputs[key_input]),
            *self._build_options(),
            f'> {path_tsv_out}'
        ])
        self._execute_command()
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]
        self._informs['cols'] = ['identity', 'hashes', 'm_mult', 'pval', 'q_id', 't']

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
