from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsCsq(BcftoolsBase):
    """
    Bcftools csq is a Haplotype-aware consequence caller.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools csq', '1.17')

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('VCF', 'VCF_GZ')):
            raise InvalidToolInputError("VCF/VCF_GZ input is required.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self._tool_outputs['VCF'] = [ToolIOFile(self._get_output_path())]

    def __build_command(self) -> None:
        """
        Builds the command line command.
        :return: None
        """
        parts = [
            self._tool_command,
            ' '.join(self._build_options()),
            str(next(self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ') if k in self._tool_inputs)),
        ]
        if 'GFF' in self._tool_inputs:
            parts.insert(3, f"--gff-annot {self._tool_inputs['GFF'][0].path}")
        if 'FASTA' in self._tool_inputs:
            parts.insert(3, f"--fasta-ref {self._tool_inputs['FASTA'][0].path}")
        self._command.command = ' '.join(parts)
