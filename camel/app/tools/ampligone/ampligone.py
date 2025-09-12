import re
from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class AmpliGone(Tool):
    """
    AmpliGone is a tool that accurately finds and removes primer sequences from NGS reads in an amplicon experiment.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('AmpliGone', '1.3.0')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA_primers' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA file with primers sequences is required")
        if 'FASTA_ref' not in self._tool_inputs:
            raise InvalidToolInputError("Reference FASTA file is required")
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError("FASTQ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Output file paths
        path_bed_out = self.folder / 'primer_locations.bed'
        fq_name_out = re.sub(r'\.fastq(\.gz)?|\.fq(\.gz)?', '_clipped.fastq', self._tool_inputs['FASTQ'][0].path.name)
        path_fq_out = self.folder / fq_name_out

        # Create & execute command
        self._build_command(path_fq_out, path_bed_out)
        self._execute_command()

        # Collect output
        self._tool_outputs['FASTQ'] = [ToolIOFile(path_fq_out)]
        self._tool_outputs['BED'] = [ToolIOFile(path_bed_out)]

        # Informs
        self._collect_informs(self._command.stdout)
        self._informs['error_rate'] = float(self._parameters['error_rate'].value)
        self._informs['fasta_primers'] = self._tool_inputs['FASTA_primers'][0].path.name

    def _build_command(self, path_fq_out: Path, path_bed_out: Path) -> None:
        """
        Creates the command to run the tool.
        :param path_fq_out: Path for trimmed/clipped FASTQ files
        :param path_bed_out: Path for BED file with primer locations
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"--reference {self._tool_inputs['FASTA_ref'][0].path}",
            f"--primers {self._tool_inputs['FASTA_primers'][0].path}",
            f"--input {self._tool_inputs['FASTQ'][0].path}",
            f"--output {path_fq_out}",
            f'--export-primers {path_bed_out}',
            *self._build_options()
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _collect_informs(self, stdout: str) -> None:
        """
        Collects the informs from the stdout.
        :param stdout: stdout
        :return: None
        """
        nucleotides_removed = None
        percentage_removed = None
        for line in stdout.splitlines():
            # Nucleotides removed
            m = re.search(r'Removed a total of (\d+) nucleotides.', line)
            if m:
                nucleotides_removed = int(m.group(1))

            # Percentage removed
            m = re.search(r'This is (\d+.\d+)% of the total amount', line)
            if m:
                percentage_removed = float(m.group(1))
        self._informs['nucleotides_removed'] = nucleotides_removed
        self._informs['percentage_removed'] = percentage_removed
