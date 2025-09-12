import json
from pathlib import Path
from typing import Any

import pandas as pd
from pandas.errors import EmptyDataError

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class Nextclade3(Tool):
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
        super().__init__('Nextclade', '3.1.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['FASTA', 'DB'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_tsv_out = self.folder / 'nextclade.tsv'
        self._command.command = ' '.join([
            self._tool_command,
            'run',
            '--input-dataset', str(self._tool_inputs['DB'][0].path),
            '--output-all', str(self.folder),
            '--output-tsv', str(path_tsv_out),
            str(self._tool_inputs['FASTA'][0].path)
        ])
        self._execute_command()
        self._parse_output_tsv(path_tsv_out)
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]
        self._informs['db'] = Nextclade3._get_database_info(self._tool_inputs['DB'][0].path)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, self._command, exit_code=0)

    def _parse_output_tsv(self, path_tsv: Path) -> None:
        """
        Parses the output TSV file and stores the output in the informs.
        :param path_tsv: Path to input TSV file
        :return: None
        """
        data_in = pd.read_table(path_tsv)
        self._informs['results'] = []
        for row in data_in.to_dict('records'):
            self._informs['results'].append(row)

    @staticmethod
    def _get_database_info(dir_db: Path) -> dict[str, Any]:
        """
        Returns the Nextclade database information.
        :return: Database information
        """
        # Parse metadata file
        path_info = dir_db / 'pathogen.json'
        logger.info(f'Retrieving DB information from: {path_info}')
        with path_info.open() as handle:
            data = json.load(handle)

        # Parse metadata columns
        path_metadata = dir_db / 'report_clade_cols.tsv'
        logger.info(f'Parsing metadata columns from: {path_metadata}')
        try:
            data_db = pd.read_table(dir_db / 'report_clade_cols.tsv')
        except EmptyDataError:
            logger.info('No metadata columns found')
            data_db = None

        # Return information
        return {
            'version': str(data['version']['tag']),
            'reference': str(data['attributes']['reference name']),
            'metadata_columns': data_db.to_dict('records') if data_db is not None else None
        }
