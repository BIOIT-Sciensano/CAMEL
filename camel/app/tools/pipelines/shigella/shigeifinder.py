import pandas as pd

from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.tool import Tool
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.errors import ToolExecutionError


class ShigEiFinder(Tool):
    """
    ShigEiFinder differentiates between Shigella/EIEC using cluster-specific genes and
    identifies the serotype using O-antigen/H-antigen genes.
    """

    def __init__(self) -> None:
        """
        Initializes the ShigEiFinder tool.
        """
        super().__init__('ShigEiFinder', '1.3.5')

    def _check_input(self) -> None:
        """
        Checks whether the provided input is valid:
        - FASTA is the only required input
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError('Only a single FASTA file can be analyzed at a time')
        super()._check_input()

    def __build_command(self, input_fasta: Path, output_tsv: Path) -> None:
        """
        Concatenates required parameters and options to build the command.
        :param input_fasta: Path to assembly file
        :param output_tsv: Path of output file
        :return: None
        """
        self._command.command = ' '.join(
            [
                self._tool_command,
                f'-i {input_fasta}',
                f'--output {output_tsv}',
                *self._build_options(),
            ]
        )

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Symlink the input FASTA file
        fasta_in = self._tool_inputs['FASTA'][0].path

        # Prepare output
        tsv_out = self.folder / 'shigeifinder_out.tsv'

        # Run the command
        self.__build_command(fasta_in, tsv_out)
        self._execute_command()

        # Collect the output
        if not tsv_out.exists():
            raise ToolExecutionError(self.name, f'{tsv_out} not generated (TSV)')
        self._tool_outputs['TSV'] = [ToolIOFile(tsv_out)]

        # Parse TSV output file
        self._parse_tsv(self._tool_outputs['TSV'][0].path)

    def __extract_species(self, serotype_abbrev: str) -> str:
        """
        Parse the output and extracts the isolate species.
        :param serotype_abbrev: Serotype abbreviation output by ShigEiFinder
        :return: Shigella species or E. coli serotype
        """
        if serotype_abbrev.startswith('SB'):
            return 'Shigella boydii'

        if serotype_abbrev.startswith('SD'):
            return 'Shigella dysenteriae'

        if serotype_abbrev.startswith('SF'):
            return 'Shigella flexneri'

        if serotype_abbrev.startswith('SS'):
            return 'Shigella sonnei'

        if serotype_abbrev.startswith('EIEC'):
            coli_serotype = serotype_abbrev.split()[1]
            return f'Enteroinvasive Escherichia coli {coli_serotype}'

        if serotype_abbrev.startswith('Not Shigella/EIEC'):
            return 'Not Shigella/EIEC'

        else:
            return str(serotype_abbrev)

    def _parse_tsv(self, path_tsv: Path) -> None:
        """
        Parses the output TSV file and stores the results in the informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data_serotype = pd.read_table(path_tsv)
        output_dict = data_serotype.fillna('-').to_dict('records')[0]
        self._informs['species'] = self.__extract_species(output_dict['SEROTYPE'])
        self._informs['serotype'] = output_dict['SEROTYPE']
        self._informs['O_antigen'] = output_dict['O_ANTIGEN']
        self._informs['H_antigen'] = output_dict['H_ANTIGEN']
        self._informs['cluster'] = output_dict['CLUSTER']
