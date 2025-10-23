from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.gatk.gatk import GATK


class GATKAnalyzeCovariates(GATK):
    """
    ==========================
    GATK AnalyzeCovariates 3.7
    ==========================
    Analyses and generates plots about covariates analysed during Base Quality Score Recalibration step and used to recalibrate bases.
    web: https://software.broadinstitute.org/gatk/gatkdocs/3.5-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php

    Required inputs:
    ----------------
    'TXT_TABLE_BEFORE': ToolIOFile object. Base quality score recalibration table on original data. GATK text format.
    'TXT_TABLE_AFTER':  ToolIOFile object. Base quality score recalibration table on recalibrated data. GATK text format.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Optional input:
    ---------------
    'BQSR':             ToolIOFile object. Base quality score recalibration table on extra data. GATK text format.

    Output:
    -------
    'PDF':              ToolIOFile object. pdf document with graphs of recalibration.
    'CSV':              ToolIOFile object. csv with recalibration covariates data used for generating the pdf output.

    Mandatory parameters:
    ---------------------
    - pdf_output       pdf output file name. Default value: 'recal_QC_plots.pdf'
    """

    def __init__(self):
        """
        Initialize GATKAnalyzeCovariates tool.
        :return: None
        """
        super().__init__('gatk AnalyzeCovariates', '3.7')

        self._required_inputs = ['TXT_TABLE_BEFORE', 'TXT_TABLE_AFTER', 'FASTA_REF']

    def _set_input(self):
        """
        Set the input specification in the input_string:
        - before and after covariates tables (and optional BQSR table)
        - reference fasta (default or superseded).
        Overrides base class method.
        :return: None
        """
        self._input_string += "-before {} ".format(self._tool_inputs['TXT_TABLE_BEFORE'][0].path)
        self._input_string += "-after {} ".format(self._tool_inputs['TXT_TABLE_AFTER'][0].path)

        self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

        if 'BQSR' in self._tool_inputs:
            self._input_string += "-BQSR {} ".format(self._tool_inputs['BQSR'][0].path)

    def _set_output(self):
        """
        Set the output specifications in the ouptut_string:
        - plots pdf file,
        - data csv file (optional)
        Supersedes _set_output in GATK class.
        :return: None
        """
        self._tool_outputs['PDF'] = [ToolIOFile(self._folder / self.get_param_value('pdf_output'))]
        self._tool_outputs['CSV'] = [ToolIOFile(self._folder / self.get_param_value('csv_output'))]
