from camel.app.error import InvalidToolInputError
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4BaseRecalibrator(GATK4):
    """
    =============================
    GATK BaseRecalibrator 4.1.9.0
    =============================
    Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.

    Required inputs:
    ----------------
    'BAM':              ToolIOFile object. Input BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Optional input:
    ---------------
    'VCF_KNOWN_SNPS':   ToolIOFile object. GATK high confidence SNP vcf file location.
    'VCF_KNOWN_INDELS': ToolIOFile object. GATK high confidence indels vcf file location.
    'TXT_intervals:     ToolIOFile object. GATK-style intervals file, BED file or Picard-style intervals list

    Output:
    -------
    'TXT_RecalibrationTable': ToolIOFile object. Text file containing recalibration data.

    Mandatory parameters:
    ---------------------
    - recal_table_output       recalibration table name. Default value: 'recalibrationData.tabl'
    """

    def __init__(self) -> None:
        """
        Initialize GATKBaseRecalibrator tool.
        :return: None
        """
        super().__init__('gatk4 BaseRecalibrator', '4.1.9.0')
        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'TXT_RecalibrationTable'

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

        if 'VCF_KNOWN_SNPS' in self._tool_inputs:
            for file in self._tool_inputs['VCF_KNOWN_SNPS']:
                self._input_string += f"--known-sites {file.path} "

        if 'VCF_KNOWN_INDELS' in self._tool_inputs:
            for file in self._tool_inputs['VCF_KNOWN_INDELS']:
                self._input_string += f"--known-sites {file.path} "

        if 'TXT_intervals' in self._tool_inputs:
            for interval in self._tool_inputs['TXT_intervals']:
                self._input_string += f"--intervals {interval.path} "

    def _check_input(self) -> None:
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        if not ('VCF_KNOWN_SNPS' in self._tool_inputs or 'VCF_KNOWN_INDELS' in self._tool_inputs):
            raise InvalidToolInputError(f'GATK {self._name} required VCF_KNOWN_SNPS or VCF_KNOWN_INDELS input '
                                                 f'is missing in _tool_inputs! (at least one should be specified)')
        super()._check_input()
