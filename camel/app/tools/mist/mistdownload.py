from camel.app.core.command import Command
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils


class MiSTDownload(Tool):
    """
    MiST is a rapid, accurate and flexible (core-genome) multi-locus sequence typing (MLST) allele caller.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MiST download', None)

    def get_version(self) -> str:
        """
        Returns the tool version.
        :return: Tool version
        """
        command = Command('mist --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        path_out = self._folder / self.get_param_value('output')
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options()
        ])
        self._execute_command()
        self._tool_outputs['DIR'] = [ToolIODirectory(path_out)]
        self._tool_outputs['TXT'] = [ToolIOFile(path_out / 'fasta_list.txt')]
        if self.get_param_value('include_profiles') is True:
            self._tool_outputs['TSV'] = [ToolIOFile(path_out / 'profiles.tsv')]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
