import json
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.tool import Tool


class MiSTLincode(Tool):
    """
    MiST is a rapid, accurate and flexible (core-genome) multi-locus sequence typing (MLST) allele caller.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MiST lincode', None)

    def get_version(self) -> str:
        """
        Returns the tool version.
        :return: Tool version
        """
        command = Command('mist --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['JSON', 'DB'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        path_out = self._folder / self.get_param_value('output')
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['JSON'][0].path),
            '--db', str(self._tool_inputs['DB'][0].path),
            *self._build_options()
        ])
        self._execute_command()
        self._tool_outputs['JSON'] = [ToolIOFile(path_out)]
        self._parse_output(path_out)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _parse_output(self, path_json: Path) -> None:
        """
        Parses the JSON output file and adds statistics to the output.
        :param path_json: Path to the JSON output file
        :return: None
        """
        with path_json.open() as handle:
            data = json.load(handle)
        required_keys = ['fields', 'lincode_full', 'lincode_partial', 'st', 'nb_loci', 'nb_matches']
        for key in required_keys:
            if not key in data:
                raise ValueError(f"The output file does not contain the field: {key}")
        self.informs.update(data)
