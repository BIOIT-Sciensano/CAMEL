from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4AnalyzeCovariates(GATK4):
    """
    ==============================
    GATK AnalyzeCovariates 4.1.9.0
    ==============================
    Analyses and generates plots about covariates analysed during Base Quality Score Recalibration step and used to recalibrate bases.

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
    - plots-report-file       pdf output file name. Default value: 'recal_QC_plots.pdf'
    """

    def __init__(self) -> None:
        """
        Initialize GATKAnalyzeCovariates tool.
        :return: None
        """
        super().__init__('gatk4 AnalyzeCovariates', '4.1.9.0')

        self._required_inputs = ['TXT_TABLE_BEFORE', 'TXT_TABLE_AFTER']

    def _set_input(self) -> None:
        """
        Set the input specification in the input_string:
        - before and after covariates tables (and optional BQSR table)
        - reference fasta (default or superseded).
        Overrides base class method.
        :return: None
        """
        self._input_string += f"--before-report-file {self._tool_inputs['TXT_TABLE_BEFORE'][0].path} "
        self._input_string += f"--after-report-file {self._tool_inputs['TXT_TABLE_AFTER'][0].path} "

        if 'BQSR' in self._tool_inputs:
            self._input_string += f"--bqsr-recal-file {self._tool_inputs['BQSR'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specifications in the ouptut_string:
        - plots pdf file,
        - data csv file (optional)
        Supersedes _set_output in GATK class.
        :return: None
        """
        self._tool_outputs['PDF'] = [ToolIOFile(Path(self._folder) / self._parameters['pdf_output'].value)]
        self._tool_outputs['CSV'] = [ToolIOFile(Path(self._folder) / self._parameters['csv_output'].value)]
