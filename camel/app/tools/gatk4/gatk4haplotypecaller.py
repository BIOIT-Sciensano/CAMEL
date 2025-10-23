from camel.app.core.errors import InvalidParameterError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4HaplotypeCaller(GATK4):
    """
    Class for GATK4 HaplotypeCaller function
    Call germline SNPs and indels via local re-assembly of haplotypes

    Required inputs:
    ----------------
    'BAM':              ToolIOFile object. Input BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Output:
    -------
    'VCF':              ToolIOFile object. (g)VCF file

    """

    def __init__(self) -> None:
        """
        Initialize the GATK4 HaplotypeCaller
        :return: None
        """
        super().__init__('gatk4 HaplotypeCaller', '4.1.9.0')
        self._required_inputs = ['BAM', 'FASTA_REF']
        self._specific_parameters = ['gqb', 'annotation_group']
        self._output_type = 'VCF'

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super()._set_input()
        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += f" --input {bam_file}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        super()._set_output()
        if 'bam_output' in self._parameters:
            self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self._parameters['bam_output'].value)]

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'gqb' in self._parameters:
            gqb_list = self._parameters['gqb'].value.split(",")
            if sorted(gqb_list) != gqb_list:
                raise InvalidParameterError("GQ Bands list must be specified in increasing order.")
        super()._check_input()

    def _build_command(self) -> None:
        """
        Build the command to run tool. Supersedes that of parent class.
        :return: None
        """
        if 'gqb' in self._parameters:
            self._option_string += self.__split_multi_options('gqb')

        if 'annotation_group' in self._parameters:
            self._option_string += self.__split_multi_options('annotation_group')
        super()._build_command()

    def __split_multi_options(self, option) -> str:
        """
        Multiple values allowed for certain parameters. These are passed in a comma separated string and need to be split
        """
        option_list = self._parameters[option].value.split(",")
        return "".join(f" {self._parameters[option].option} {s} " for s in option_list)
