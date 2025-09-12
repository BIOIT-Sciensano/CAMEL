from pathlib import Path

import pandas as pd

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Nextclade(Tool):
    """
    Nextclade is a tool that identifies differences between your sequences and a reference sequence, uses these
    differences to assign your sequences to clades, and reports potential sequence quality issues in your data. You can
    use the tool to analyze sequences before you upload them to a database, or if you want to assign Nextstrain clades
    to a set of sequences.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
                :return: None
        """
        super().__init__('Nextclade', '2.14.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('Database input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            'run',
            '--input-dataset', str(self._tool_inputs['DB'][0].path),
            '--output-all', str(self.folder),
            str(self._tool_inputs['FASTA'][0].path)
        ])
        self._execute_command()
        path_csv_out = self.folder / 'nextclade.csv'
        self._parse_output_csv(path_csv_out)
        self._tool_outputs['CSV'] = [ToolIOFile(path_csv_out)]

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _parse_output_csv(self, path_csv: Path) -> None:
        """
        Parses the output CSV file and stores the output in the informs.
        :param path_csv: Path to input CSV file
        :return: None
        """
        data_in = pd.read_table(path_csv, sep=';')
        self._informs['results'] = []
        for row in data_in.to_dict('records'):
            self._informs['results'].append(row)
