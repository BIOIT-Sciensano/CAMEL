import json
from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class MiSTCall(Tool):
    """
    MiST is a rapid, accurate and flexible (core-genome) multi-locus sequence typing (MLST) allele caller.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MiST call', None)

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
        toolutils.check_input(self, keys_required=['DB', 'FASTA'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        path_out = self._folder / 'mist_out.json'
        self._command.command = ' '.join([
            self._tool_command,
            '--fasta', str(self._tool_inputs['FASTA'][0].path),
            '--db', str(self._tool_inputs['DB'][0].path),
            '--out-json', str(path_out),
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
        loci_detected = sum(res['allele_str'] != '-' for _, res in data['alleles'].items())
        self._informs['loci_detected'] = loci_detected
        self._informs['loci_total'] = len(data['alleles'])
