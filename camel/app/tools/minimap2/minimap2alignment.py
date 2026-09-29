from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.minimap2.minimap2base import Minimap2Base


class Minimap2Alignment(Minimap2Base):
    """
    Aligns query sequences against target sequences and reports the alignments in PAF format.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Minimap2 alignment', version=None)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['FASTA', 'FASTA_REF'], keys_allowed=['FASTA', 'FASTA_REF'])
        if len(self._tool_inputs['FASTA_REF']) != 1:
            raise InvalidToolInputError('Only one FASTA_REF input file is supported')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / self.get_param_value('output_filename')
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(excluded_parameters=['output_filename']),
            str(self._tool_inputs['FASTA_REF'][0].path),
            *(str(x.path) for x in self._tool_inputs['FASTA']),
            f'> {path_out}'
        ])
        self._execute_command()
        self._tool_outputs['PAF'] = [ToolIOFile(path_out)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
