from pathlib import Path

import pandas as pd
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool


class Pharokka(Tool):
    """
    Pharokka is a fast phage genome annotation tool.
    """

    OUTPUT_DICT = {
        'GBK': 'pharokka.gbk',
        'TSV_CARD': 'top_hits_card.tsv',
        'TSV_VFDB': 'top_hits_vfdb.tsv',
        'TSV_STATS': 'pharokka_cds_functions.tsv',
        'TSV_INPHARED': 'pharokka_top_hits_mash_inphared.tsv'
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Pharokka', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.strip()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()

        # Collect the output
        for key, basename in Pharokka.OUTPUT_DICT.items():
            path_out = Path(self.folder, self._parameters['outdir'].value, basename)
            if not path_out.exists():
                raise ToolExecutionError(self.name, f'{path_out} not generated ({key})')
            self._tool_outputs[key] = [ToolIOFile(path_out)]

        # Parse output files
        self._parse_card_hits(self._tool_outputs['TSV_CARD'][0].path)
        self._parse_vfdb_hits(self._tool_outputs['TSV_VFDB'][0].path)
        self._parse_metrics(self._tool_outputs['TSV_STATS'][0].path)

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("No FASTA input found")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError("Only one FASTA file can be annotated at a time.")
        super()._check_input()

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __build_command(self) -> None:
        """
        Build the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            '-i',
            str(self._tool_inputs['FASTA'][0].path)
        ])

    def _parse_metrics(self, path_tsv: Path) -> None:
        """
        Parses the file with the summary metrics and stores them in informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data = pd.read_table(path_tsv)
        cds_count = data.iloc[0]['Count']
        self._informs['stats'] = f'{cds_count} CDS were annotated by Pharokka.'

    def _parse_card_hits(self, path_tsv: Path) -> None:
        """
        Parses the top CARD hits TSV file and stores the results in the informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data = pd.read_table(path_tsv)
        if data.empty:
            self._informs['card_hits'] = 'No AMR genes were detected.'
        else:
            self._informs['card_hits'] = f'{len(data)} AMR genes were detected.'

    def _parse_vfdb_hits(self, path_tsv: Path) -> None:
        """
        Parses the top VFDB hits TSV file and stores the results in the informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data = pd.read_table(path_tsv)
        if data.empty:
            self._informs['vfdb_hits'] = 'No virulence genes were detected.'
        else:
            self._informs['vfdb_hits'] = f'{len(data)} virulence genes were detected.'
