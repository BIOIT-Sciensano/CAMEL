import time
from pathlib import Path

import pandas as pd
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class CharacterizeNeisseriaCapsule(Tool):
    """
    characterize_neisseria_capsule is a tool implementing a WGS-based method for N. meningitidis
    serogroup predictions by identifying capsule genes and genetic variations that might impact their expression.
    """

    MAX_ATTEMPTS = 3
    RETRY_DELAY = 3
    PUBMLST_ERRORS = ('JSONDecodeError', 'MaxRetryError')

    def __init__(self) -> None:
        """
        Initializes the characterize_neisseria_capsule tool.
        """
        super().__init__('characterize_neisseria_capsule', 'a75a009')

    def _check_input(self) -> None:
        """
        Checks whether the provided inputs is valid:
        - FASTA is the only required input
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Only a single FASTA file can be analyzed at a time')
        super()._check_input()

    def __build_command(self, dir_path: Path, dir_out: Path) -> None:
        """
        Concatenates required parameters and options to build the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f'--out {dir_out}',
            f'--indir {dir_path}',
            *self._build_options()
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __is_pubmlst_error(self) -> bool:
        """
        Checks if the last execution failed on a PubMLST query (invalid JSON response or connection error).
        :return: True if the failure is caused by PubMLST, False otherwise
        """
        stderr = self._command.stderr or ''
        return any(error in stderr for error in CharacterizeNeisseriaCapsule.PUBMLST_ERRORS)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Symlink the input FASTA file
        dir_fasta_in = self.folder / 'fasta_in'
        dir_fasta_in.mkdir()
        (dir_fasta_in / self._tool_inputs['FASTA'][0].path.name).symlink_to(self._tool_inputs['FASTA'][0].path)

        # Run the command, retrying when the PubMLST allele lookup fails (e.g. rate limiting). The script reuses the
        # BLAST results in the existing output folder, so a retry only repeats the PubMLST queries.
        dir_out = self.folder / 'out'
        for attempt in range(1, CharacterizeNeisseriaCapsule.MAX_ATTEMPTS + 1):
            self.__build_command(dir_fasta_in, dir_out)
            try:
                self._execute_command()
                break
            except ToolExecutionError:
                if attempt == CharacterizeNeisseriaCapsule.MAX_ATTEMPTS or not self.__is_pubmlst_error():
                    raise
                delay = CharacterizeNeisseriaCapsule.RETRY_DELAY * 2 ** (attempt - 1)
                logger.warning(
                    f'PubMLST query failed (attempt {attempt}/{CharacterizeNeisseriaCapsule.MAX_ATTEMPTS}), '
                    f'retrying in {delay} seconds'
                )
                time.sleep(delay)

        # Collect the output
        try:
            self._tool_outputs['TSV'] = [ToolIOFile(next((dir_out / 'serogroup').glob('serogroup_predictions_*.tab')))]
        except StopIteration:
            raise ToolExecutionError(self.name, f"TSV file not found in output folder: {dir_out / 'serogroup'}")
        self._tool_outputs['JSON'] = [ToolIOFile(dir_out / 'serogroup' / 'serogroup_results.json')]
        self._parse_tsv(self._tool_outputs['TSV'][0].path)

    def _parse_tsv(self, path_tsv: Path) -> None:
        """
        Parses the output TSV file and stores the results in the informs.
        :return: None
        """
        data_sero = pd.read_table(path_tsv)
        output_dict = data_sero.where(data_sero.notna(), '-').to_dict('records')[0]
        self._informs['detected_serogroup'] = output_dict['SG']
        self._informs['genes_present'] = output_dict['Genes_Present']
        self._informs['notes'] = output_dict['Notes']
