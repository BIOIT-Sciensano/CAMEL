from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.gatk.gatk import GATK


class GATKBaseRecalibrator(GATK):
    """
    ==========================
    GATK BaseRecalibrator 3.7
    ==========================
    Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.
    web: https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
    https://software.broadinstitute.org/gatk/gatkdocs/4.beta.5/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
    https://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

    Required inputs:
    ----------------
    'BAM':              ToolIOFile object. Input BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Optional input:
    ---------------
    'VCF_KNOWN_SNPS':   ToolIOFile object. GATK high confidence SNP vcf file location.
    'VCF_KNOWN_INDELS': ToolIOFile object. GATK high confidence indels vcf file location.

    Output:
    -------
    'TXT_RecalibrationTable': ToolIOFile object. Text file containing recalibration data.

    Mandatory parameters:
    ---------------------
    - recal_table_output       recalibration table name. Default value: 'recalibrationData.tabl'
    """

    def __init__(self):
        """
        Initialize GATKBaseRecalibrator tool.
        :return: None
        """
        super().__init__('gatk BaseRecalibrator', '3.7')

        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'TXT_RecalibrationTable'
        self.__snps_known_sites_path = ''
        self.__indels_known_sites_path = ''

    def _set_input(self):
        """
        Set the input specification in the input_string
        Overrides method in parent class.
        :return: None
        """
        # set input BAM
        self._input_string += "-I {} ".format(self._tool_inputs['BAM'][0].path)

        # set reference genome, known snps and known indels
        self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

        if 'VCF_KNOWN_SNPS' in self._tool_inputs:
            self._input_string += "-knownSites {} ".format(self._tool_inputs['VCF_KNOWN_SNPS'][0].path)

        if 'VCF_KNOWN_INDELS' in self._tool_inputs:
            self._input_string += "-knownSites {} ".format(self._tool_inputs['VCF_KNOWN_INDELS'][0].path)

        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

    def _set_output(self):
        """
        Set the output specification in the output_string
        Overrides method in parent class.
        :return: None
        """
        self._tool_outputs['TXT_RecalibrationTable'] = [
            ToolIOFile(self._folder / self.get_param_value('recal_table_output'))]
