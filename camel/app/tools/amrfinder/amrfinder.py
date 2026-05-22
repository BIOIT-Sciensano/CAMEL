from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class AMRFinder(Tool):
    """
    NCBI Antimicrobial Resistance Gene Finder (AMRFinderPlus).
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('AMRFinder', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd = True)
        return command.stdout.strip()

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DIR)")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = Path(self.folder) / Path(self._parameters['output_path'].value)
        dir_temp = config.dir_temp
        self._command.command = ' '.join([
            self._tool_command,
            '--nucleotide', str(self._tool_inputs['FASTA'][0].path),
            '--database', str(self._tool_inputs['DIR'][0].path),
            str(self._parameters['output_path'].option), str(output_path)
        ] + self._build_options(['output_path']))
        self._execute_command(env={'TMPDIR': str(dir_temp)})
        self._tool_outputs['TSV'] = [ToolIOFile(output_path)]
        self._informs['db_version'] = self._tool_inputs['DIR'][0].path.resolve().name

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
