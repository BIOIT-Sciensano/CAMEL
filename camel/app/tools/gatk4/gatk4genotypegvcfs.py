from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4GenotypeGVCFs(GATK4):
    """
    Class for GATK GenotypeGVCFs function
    """

    def __init__(self) -> None:
        """
        Initialize the GATK4GenotypeGVCFs tool
        :return: None
        """
        super().__init__('gatk4 GenotypeGVCFs', '4.1.9.0')

        self._required_inputs = ['gVCF', 'FASTA_REF']
        self._output_type = 'VCF_MultipleSample'
        self._specific_parameters = ['gendb']

    def _check_input(self) -> None:
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        if self._parameters['gendb'] is not None:
            self._required_inputs.remove('gVCF')
        super()._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super()._set_input()
        if 'gVCF' in self._tool_inputs:
            for f in self._tool_inputs['gVCF']:
                self._input_string += f"--variant {f.path} "
        elif 'gendb' in self._parameters:
            self._input_string += f"--variant gendb:/{self._parameters['gendb'].value}"

    def _build_command(self):
        """
        Build the command to run tool.
        :return: None
        """
        if 'annotation_group' in self._parameters:
            self._option_string += ' --annotation-group '
            self._option_string += ' --annotation-group '.join(self._parameters['annotation_group'].value.split(","))
        self._option_string += ' '
        self._specific_parameters.append("annotation_group")
        super()._build_command()
