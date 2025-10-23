from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class ALE2Wiggle(Tool):
    """
    ALE is the Assembly Likelihood Evaluation framework that systematically evaluates the accuracy of an assembly
    in a reference-independent manner using rigorous statistical methods.

    Ale2Wiggle converts the ALE output file to a set of wiggle files, which can be opened in IGV.
    """

    def __init__(self) -> None:
        """
        Initialize the ALE2Wiggle tool.
        :return: None
        """
        super().__init__('ALE2Wiggle', '2022.05.03')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        try:
            path_symlink = self._folder / Path(str(self._tool_inputs['ALE'][0])).name
            path_symlink.symlink_to(self._tool_inputs['ALE'][0].path)
            self._build_command(path_symlink)
        except FileExistsError:
            self._build_command(self._tool_inputs['ALE'][0].path)

        self._execute_command()
        self._set_output()

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'ALE' not in self._tool_inputs:
            raise InvalidToolInputError('ALE file is required')
        super()._check_input()

    def _build_command(self, ale_output: Path) -> None:
        """
        Builds the command to run ALE2Wiggle.
        :ale_output: Path to the file output by ALE
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, str(ale_output)])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks command output.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool outputs.
        """
        basename = Path(str(self._tool_inputs['ALE'][0])).name
        tsv_outputs = [f'{basename}-{stat}.wig' for stat in ['depth', 'kmer', 'insert', 'place']]
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / f'{f}') for f in tsv_outputs]
