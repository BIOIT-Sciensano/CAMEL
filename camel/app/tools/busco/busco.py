import json
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Busco(Tool):
    """
    Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs, BUSCO metric
    is complementary to technical metrics like N50.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('busco', '5.5.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('DB input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_out = Path(self._folder) / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            f"--in {self._tool_inputs['FASTA'][0].path}",
            f'--out {dir_out.name}',
            f"--download_path {self._tool_inputs['DB'][0].path}",
            '--mode', 'genome',
            *self._build_options(),
        ])
        self._execute_command()
        self._tool_outputs['TXT'] = [ToolIOFile(next(p for p in dir_out.glob('short_summary.*.txt')))]
        self._set_informs(dir_out)

    def _set_informs(self, dir_out: Path) -> None:
        """
        Collects the informs from the output JSON file.
        :param: dir_out: Output directory
        :return: None
        """
        path_json = next(dir_out.glob('short_summary.*.json'))
        with path_json.open() as handle:
            self._informs['results'] = json.load(handle)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
