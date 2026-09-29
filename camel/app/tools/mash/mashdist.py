from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.mash.mashbase import MashBase


class MashDist(MashBase):
    """
    Mash dist estimates the distance of each query sequence to each reference sequence.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('mash dist', version=None)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['DB', 'FASTA'], keys_allowed=['DB', 'FASTA'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_tsv_out = self.folder / 'mash_dist.tsv'
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            str(self._tool_inputs['DB'][0].path),
            *(str(x.path) for x in self._tool_inputs['FASTA']),
            f'> {path_tsv_out}'
        ])
        self._execute_command()
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]
        self._informs['cols'] = ['reference_id', 'query_id', 'distance', 'pval', 'shared_hashes']

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
