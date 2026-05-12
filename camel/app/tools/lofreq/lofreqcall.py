from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.lofreq.lofreq import Lofreq


class LofreqCall(Lofreq):
    """
    LoFreq is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data.
    """

    def __init__(self) -> None:
        """
        Initializes Lofreq call.
        :return: None
        """
        super().__init__('Lofreq call', version=None)

    def _check_input(self) -> None:
        """
        Checks that the input is correct.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA reference is required')
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError('BAM alignment file is required')
        if len(self._tool_inputs['BAM']) != 1:
            raise ValueError("Exactly one BAM input file expected")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes Lofreq call.
        :return: None
        """
        fasta_input = Path(str(self._tool_inputs['FASTA'][0]))
        bam_input = Path(str(self._tool_inputs['BAM'][0]))
        self.__build_command(fasta_input, bam_input)
        self._execute_command()
        self.__set_output()

    def __build_command(self, fasta_input: Path, bam_input: Path) -> None:
        """
        Builds the command.
        :param fasta_input: Path to the FASTA input file
        :param bam_input: Path to the BAM input file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, f'--ref {fasta_input}',
            *self._build_options(),
            str(bam_input)
        ])

    def __set_output(self) -> None:
        """
        Sets the output of Lofreq call.
        :return: None
        """
        output_file_path = self.folder / self._parameters['output_filename'].value
        self._tool_outputs['VCF'] = [ToolIOFile(output_file_path)]
