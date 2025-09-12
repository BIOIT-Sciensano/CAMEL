from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GrapeTree(Tool):
    """
    Package to create cgMLST-based phylogenies associated with EnteroBase.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('GrapeTree', '2.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("Allele matrix (TSV) input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self.folder / self._parameters['output_path'].value
        self._command.command = ' '.join([
            self._tool_command,
            '--profile', str(self._tool_inputs['TSV'][0].path),
            *self._build_options(excluded_parameters=['output_path']),
            f'> {output_path}'
        ])
        self._execute_command()
        self._tool_outputs['NWK'] = [ToolIOFile(output_path)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
