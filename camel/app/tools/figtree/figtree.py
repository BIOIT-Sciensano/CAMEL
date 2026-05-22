from pathlib import Path

from Bio import Phylo
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class FigTree(Tool):
    """
    FigTree is designed as a graphical viewer of phylogenetic trees.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('FigTree', '1.4.4')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'NWK' not in self._tool_inputs:
            raise InvalidToolInputError('Newick input is required (NWK)')
        super()._check_input()

    def _get_tool_command(self) -> str:
        """
        Returns the tool command.
        :return: String representing the command
        """
        if config.dependency_service == 'lmod':
            return 'java -jar $JAR_FIGTREE'
        return self._tool_command

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if 'TXT' not in self._tool_inputs:
            # Add basic Newick file when no FigTree template is provided
            input_file = self._tool_inputs['NWK'][0].path
        else:
            # Add the Figtree template information to the input tree (converted to Nexus format)
            input_file = self._add_template_info(self._tool_inputs['TXT'][0].path)
        output_path = self.folder / self._parameters['output_path'].value

        # Construct and execute command
        self._command.command = ' '.join([
            self._get_tool_command(),
            '-graphic PNG',
            *self._build_options(excluded_parameters=['output_path']),
            str(input_file),
            str(output_path)
        ])
        self._execute_command()

        # Set output
        self.tool_outputs['PNG'] = [ToolIOFile(output_path)]

    def _add_template_info(self, path_template: Path) -> Path:
        """
        Adds the template information.
        :param path_template: Path to the template file
        :return: Updated input file in Nexus format
        """
        path_nexus = self.folder / 'input_tree.nexus'
        Phylo.convert(str(self._tool_inputs['NWK'][0].path), 'newick', str(path_nexus), 'nexus')
        with path_nexus.open('a', encoding='utf-8') as handle_out, path_template.open(encoding='utf-8') as handle_in:
            # Add template content
            handle_out.write(handle_in.read())
            handle_out.write('\n')
        return path_nexus

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
