from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4ApplyBQSR(GATK4):
    """
    =============================
    GATK ApplyBQSR 4.1.9.0
    =============================
    Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.

    Required inputs:
    ----------------
    'BAM':              ToolIOFile object. Input BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.
    'BQSR':             ToolIOFile object. Input recalibration table for BQSR.

    Optional inputs:
    ----------------
    'TXT_intervals'     ToolIOFile object. GATK-style intervals file, BED file or Picard-style intervals list

    Output:
    -------
    'BAM':              ToolIOFile object. Output BAM file containing recalibrated read data.

    """

    def __init__(self) -> None:
        """
        Initialize GATK4ApplyBQSR tool.
        :return: None
        """
        super().__init__('gatk4 ApplyBQSR', '4.1.9.0')
        self._required_inputs = ['BAM', 'FASTA_REF', 'BQSR']
        self._output_type = 'BAM'
        self._specific_parameters = 'static_quantized_quals_multi'

    def _set_input(self) -> None:
        """
        Set the input specification in the input_string
        Overrides method in parent class.
        :return: None
        """
        # set input BAM
        self._input_string += f"--input {self._tool_inputs['BAM'][0].path} "

        # set reference genome, known snps and known indels
        self._input_string += f"--reference {self._tool_inputs['FASTA_REF'][0].path} "

        # set input BQSR
        self._input_string += f"--bqsr-recal-file {self._tool_inputs['BQSR'][0].path} "

        # if intervals is determined, set input intervals
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += f"--intervals {self._tool_inputs['TXT_intervals'][0].path} "

    def _build_command(self) -> None:
        """
        Build the command to run tool.
        Overrides that of parent class.
        :return: None
        """
        if 'static_quantized_quals_multi' in self._parameters:
            self._option_string += ' --static-quantized-quals '
            self._option_string += ' --static-quantized-quals '.join(
                self._parameters['static_quantized_quals_multi'].value.split(","))
            self._option_string += " "
        super()._build_command()
