import re
from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.samtools.samtoolsbase import SamtoolsBase


class SamtoolsAmpliconClip(SamtoolsBase):
    """
    Clip oligos from the end of reads.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('samtools ampliconclip')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM', 'BED'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / self.get_param_value('output')
        self._command.command = ' '.join([
            self._tool_command,
            f"-b {self._tool_inputs['BED'][0].path}",
            str(self._tool_inputs['BAM'][0].path),
            *self._build_options(),
        ])
        self._execute_command()
        self._set_output(path_out)
        self._parse_stderr(self._command)

    def _set_output(self, path_out: Path) -> None:
        """
        Collects the tool output.
        :param path_out: Path to the output file
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(path_out)]

    def _parse_stderr(self, command: Command) -> None:
        """
        Parses the stderr and saves the information in the informs.
        :param command: Executed command
        :return: None
        """
        self._informs['stats'] = {}
        for line in command.stderr.splitlines():
            m = re.match(r'^([A-Z ]+): (\d+)', line)
            if not m:
                continue
            self._informs['stats'][m.group(1)] = int(m.group(2))
