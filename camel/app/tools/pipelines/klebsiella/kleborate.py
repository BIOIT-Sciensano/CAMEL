from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Kleborate(Tool):
    """
    Kleborate is a tool to screen genome assemblies of Klebsiella pneumoniae and the Klebsiella pneumoniae species
    complex (KpSC) for:
    * MLST sequence type
    * species (e.g. K. pneumoniae, K. quasipneumoniae, K. variicola, etc.)
    * ICEKp associated virulence loci: yersiniabactin (ybt), colibactin (clb), salmochelin (iro), hypermucoidy (rmpA)
    * virulence plasmid associated loci: salmochelin (iro), aerobactin (iuc), hypermucoidy (rmpA, rmpA2)
    * antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
    * K (capsule) and O antigen (LPS) serotype prediction, via wzi alleles and Kaptive
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Kleborate', version=None)

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._tool_command} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.split(' ')[-1].strip()

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / 'kleborate.txt'
        self._command.command = ' '.join([
            self._tool_command,
            f"-a {self._tool_inputs['FASTA'][0].path}",
            f'-o {path_out}',
            *self._build_options()
        ])
        self._execute_command()
        self._parse_output_file(path_out)
        self._tool_outputs['TSV'] = [ToolIOFile(path_out)]

    def _parse_output_file(self, path_out: Path) -> None:
        """
        Parses the output file and stores the results in the informs.
        :param path_out: Output file path
        :return: None
        """
        with path_out.open() as handle:
            header = handle.readline().split('\t')
            values = handle.readline().split('\t')
            for k, v in zip(header, values):
                self._informs[k] = v

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
