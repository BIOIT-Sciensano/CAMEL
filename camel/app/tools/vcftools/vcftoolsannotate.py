from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.vcftools.vcftools import VCFtools


class VCFtoolsAnnotate(VCFtools):
    """
    Annotate VCF file, add filters or custom annotations
    Reads an input VCF from stdin and prints output VCF to stdout
    """

    def __init__(self) -> None:
        """
        Initialize vcf-annotate tool.
        :return: None
        """
        super().__init__('VCFtools vcf-annotate', '0.1.16')
        self._specific_parameters = ['output']

    def _build_command(self) -> None:
        """
        Build command of vcf-annotate function
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters)

        if 'VCF' in self._tool_inputs:
            input_command = f'cat {self._tool_inputs["VCF"][0].path}'
        elif 'VCF_GZ' in self._tool_inputs:
            input_command = f'gunzip -c {self._tool_inputs["VCF_GZ"][0].path}'
        else:
            raise InvalidToolInputError("VCFtools vcf-annotate requires a VCF or VCF_GZ input file.")

        self._command.command = '|'.join([
            input_command,
            f' {self._tool_command} {" ".join(build_options)} ',
            f' bgzip -c > {self.folder / self._parameters["output"].value}'
        ])
