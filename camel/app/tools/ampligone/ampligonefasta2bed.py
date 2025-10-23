from pathlib import Path

import pandas as pd
from Bio import SeqIO

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class AmpliGoneFasta2Bed(Tool):
    """
    AmpliGone is a tool that accurately finds and removes primer sequences from NGS reads in an amplicon experiment.
    FASTA2BED determines the position of primers in the provided reference sequence and outputs a BED file with the
    corresponding coordinates.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('AmpliGone fasta2bed', '1.3.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA_primers' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA file with primers sequences is required")
        if 'FASTA_ref' not in self._tool_inputs:
            raise InvalidToolInputError("Reference FASTA file is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Output file paths
        path_bed_out = self.folder / 'primer_locations.bed'

        # Create & execute command
        self._command.command = ' '.join([
            self._tool_command,
            f"--reference {self._tool_inputs['FASTA_ref'][0].path}",
            f"--primers {self._tool_inputs['FASTA_primers'][0].path}",
            f"--output {path_bed_out}",
            *self._build_options()
        ])
        self._execute_command()

        # Collect output
        self._tool_outputs['BED'] = [ToolIOFile(path_bed_out)]

        # Informs
        self._collect_primer_stats(self._tool_inputs['FASTA_primers'][0].path, path_bed_out)
        self._informs['primer_mismatch_rate'] = float(self._parameters['primer_mismatch_rate'].value)
        self._informs['fasta_primers'] = self._tool_inputs['FASTA_primers'][0].path.name

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _collect_primer_stats(self, path_fasta_primers: Path, path_bed_out: Path) -> None:
        """
        Collects the informs from the stdout.
        :param path_fasta_primers: FASTA file with primer sequences
        :param path_bed_out: Output BED file
        :return: None
        """
        with path_fasta_primers.open() as handle:
            self._informs['primers_in'] = [s.id for s in SeqIO.parse(handle, 'fasta')]
        data_primers = pd.read_table(path_bed_out, usecols=[3], names=['primer'])
        count_by_primer = data_primers['primer'].value_counts().to_dict()
        self._informs['primers_out'] = {p: count_by_primer.get(p, 0) for p in self._informs['primers_in']}
