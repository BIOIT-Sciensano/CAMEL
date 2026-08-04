from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastautils

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Freebayes(Tool):
    """
    Freebayes is a Bayesian genetic variant detector designed to find small polymorphisms,
    specifically SNPs, indels, MNPs (multi-nucleotide polymorphisms), and complex events
    (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Freebayes', '1.3.6')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['FASTA', 'BAM'])
        if not fastautils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidToolInputError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, bam_input: Path) -> None:
        """
        Builds the command to run freebayes.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, f'--bam {bam_input}', f'--fasta-reference {fasta_input}', *self._build_options()])

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
        self._tool_outputs['VCF'] = [ToolIOFile(Path(self._folder / self._parameters['vcf'].value))]

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
