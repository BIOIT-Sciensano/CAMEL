from camel.app.core.command import Command
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.tool import Tool
from camel.app.core.utils import toolutils


class MiSTIndex(Tool):
    """
    MiST is a rapid, accurate and flexible (core-genome) multi-locus sequence typing (MLST) allele caller.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MiST index', None)

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
        if all(x not in self._tool_inputs for x in ('FASTA', 'TXT')):
            raise InvalidToolInputError("FASTA or TXT input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        dir_out = self._folder / self.get_param_value('output')
        if 'TXT' in self._tool_inputs:
            input_opt = f'--fasta-list {self._tool_inputs["TXT"][0].path}'
        else:
            input_opt = f'--fasta {" ".join(str(io.path) for io in self._tool_inputs["FASTA"])}'

        self._command.command = ' '.join([
            self._tool_command,
            input_opt,
            f'--profiles {str(self._tool_inputs["TSV"][0].path)}' if 'TSV' in self._tool_inputs else '',
            *self._build_options()
        ])
        self._execute_command()
        self._tool_outputs['DIR'] = [ToolIODirectory(dir_out)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
