from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.components.files.fastautils import FastaUtils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Sniffles(Tool):
    """
    Sniffles is a fast structural variant caller for long-read sequencing, it accurately detects SVs on germline,
    somatic and population-level for PacBio and Oxford Nanopore read data.
    """

    def __init__(self) -> None:
        """
        Initializes Sniffles 2.0.7.
                :return: None
        """
        super().__init__('Sniffles', '2.2')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA file is required')
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError('BAM alignment file is required')

        if not FastaUtils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidToolInputError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, bam_input: Path) -> None:
        """
        Builds the command to run sniffles.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, f'--input {bam_input}', f'--reference {fasta_input}', *self._build_options()])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool output.
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(self.folder / Path('variants.vcf'))]

    def _parse_output(self, path_vcf: Path) -> None:
        """
        Parses the output vcf of sniffles and store the variants found in the informs.
        Note: this method uses manual parsing because the VCF file is not compatible with PyVCF
        :path_vcf: Path to the output VCF file
        :return: None
        """
        self._informs['variants'] = {'BND': 0, 'INS': 0, 'DEL': 0, 'DUP': 0, 'INV': 0}
        with open(path_vcf) as handle:
            for line in handle.readlines():
                if not line.startswith('#'):
                    self._informs['variants'][line.split('\t')[2].split('.')[1]] += 1

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = Path(str(self._tool_inputs['FASTA'][0]))
        bam_input = Path(str(self._tool_inputs['BAM'][0]))
        self._build_command(fasta_input, bam_input)
        self._execute_command()
        self._set_output()
        self._parse_output(self.folder / Path('variants.vcf'))
