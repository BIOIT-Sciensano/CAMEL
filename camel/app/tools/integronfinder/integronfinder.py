import pandas as pd

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class IntegronFinder(Tool):
    """
    Finds integrons in DNA sequences.
    Reference: https://github.com/gem-pasteur/Integron_Finder
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('IntegronFinder', '2.0.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['FASTA'][0].path),
            *self._build_options()
        ])
        self._execute_command()
        dir_out = next(d for d in self.folder.glob('Results_*') if d.is_dir())
        path_tsv_out = next(dir_out.glob('*.integrons'))
        try:
            data_out = pd.read_table(path_tsv_out, comment='#')
            self._informs['nb_detected'] = len(data_out)
        except pd.errors.EmptyDataError as err:
            self._informs['nb_detected'] = 0
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
