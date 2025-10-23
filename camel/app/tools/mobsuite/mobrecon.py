from pathlib import Path

import pandas as pd

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class MOBRecon(Tool):
    """
    MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('MOB-recon', '3.1.8')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError("Database input (DB) is required")
        super()._check_input()

    def _build_command(self, dir_out: Path) -> None:
        """
        Builds the command line call.
        :param dir_out: Output directory
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-i {self._tool_inputs['FASTA'][0].path}",
            '--database_directory', str(self._tool_inputs['DB'][0].path),
            f'-o {dir_out}',
            *self._build_options()
        ])

    def _collect_tool_output(self, dir_out: Path) -> None:
        """
        Collects the tool output.
        :param dir_out: Output directory
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(dir_out / 'mobtyper_results.txt')]
        self._tool_outputs['TSV_contigs'] = [ToolIOFile(dir_out / 'contig_report.txt')]
        self._tool_outputs['FASTA'] = []
        for path_fasta in sorted(dir_out.glob('plasmid*.fasta')):
            self._tool_outputs['FASTA'].append(ToolIOFile(path_fasta))

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_out = self.folder / 'out'
        self._build_command(dir_out)
        self._execute_command()
        self._parse_output_file(dir_out / 'mobtyper_results.txt')
        self._parse_contig_report(dir_out / 'contig_report.txt')
        self._collect_tool_output(dir_out)

    def _parse_output_file(self, path_out: Path) -> None:
        """
        Parses the output file and stores the results in the informs.
        :param path_out: Output file path
        :return: None
        """
        self._informs['detected_plasmids'] = []
        if path_out.exists():
            data_plasmids = pd.read_table(path_out)
            logger.info(f'{len(data_plasmids)} plasmids detected')
            for record in data_plasmids.to_dict('records'):
                self._informs['detected_plasmids'].append(record)
        else:
            logger.info('No plasmids detected, creating empty output file')
            path_out.touch()

    def _parse_contig_report(self, path_out: Path) -> None:
        """
        Parses the contig report and stores the results in the informs.
        :param path_out: Output file path
        :return: None
        """
        data_in = pd.read_table(path_out)
        self._informs['contig_report'] = {
            ctg: cluster_id if cluster_id != '-' else None for
            ctg, cluster_id in zip(data_in['contig_id'], data_in['primary_cluster_id'])}

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
