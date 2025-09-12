from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class CheckV(Tool):
    """
    Assessing the quality of metagenome-assembled viral genomes
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('CheckV', '1.0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_out = self.folder / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            'end_to_end',
            str(self._tool_inputs['FASTA'][0].path),
            str(dir_out),
            *self._build_options()
        ])
        self._execute_command()
        self.__set_output(dir_out)

    def __set_output(self, dir_out: Path) -> None:
        """
        Sets the tool output for this tool.
        :param dir_out: Output directory
        :return: None
        """
        for filename in [Path(x) for x in (
                'complete_genomes.tsv', 'completeness.tsv', 'contamination.tsv', 'quality_summary.tsv')]:
            self._tool_outputs[f'TSV_{filename.stem}'] = [ToolIOFile(dir_out / filename)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
