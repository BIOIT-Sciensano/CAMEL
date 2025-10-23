from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.gatk.gatk import GATK


class GATKPrintReads(GATK):
    """
    ==========================
    GATK PrintReads 3.7
    ==========================
    Print reads from a SAM/BAM/CRAM file.
    Example: apply recalibration table to BAM file.
    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_PrintReads.php

    Required inputs:
    ----------------
    'BAM':              ToolIOFile object. Input BAM file.
    'FASTA_REF':        ToolIOFile object. FASTA file containing the reference genome.

    Optional inputs:
    ----------------
    'BQSR':             ToolIOFile object. BQSR file from base quality recalibration step to apply to input BAM.
    "TXT_intervals":    ToolIOFile object. Intervals list to restrict breadth and accelerate analysis. Bed or GATK intervals list.

    Output:
    -------
    'BAM':              ToolIOFile object. TBAM file containing (eg.) recalibrated reads.

    Parameters:
    ----------
    - bam_output       Output BAM file name. Default value: 'recalibrated.bam'
    """

    def __init__(self):
        """
        Initialize GATKPrintReads tool.
        :return: None
        """
        super().__init__('gatk PrintReads', '3.7')

        self._required_inputs = ['BAM', 'FASTA_REF']
        self._specific_parameters = ["bam_external_output"]

    def _set_input(self):
        """
        Set the input specification in the input_string
        Overrides method in parent class
        :return: None
        """
        # set input BAM
        self._input_string += "-I {} ".format(self._tool_inputs['BAM'][0].path)

        # set input recalibration table
        if 'BQSR' in self._tool_inputs:
            self._input_string += "-BQSR {} ".format(self._tool_inputs['BQSR'][0].path)

        # Intervals
        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

        # set reference genome
        self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

    def _set_output(self):
        """
        Set the output specification in the output_string:
        BAM and BAI files
        Supersedes the _set_output fct in GATK class.
        :return: None
        """
        bam_output_file = self._parameters['bam_output'].value
        self._tool_outputs['BAM'] = [ToolIOFile(self._folder / bam_output_file)]

    def run(self, folder='.'):
        """
        Runs this tool.
        :param folder: Folder to run the tool in.
        :return: None
        """
        super(GATKPrintReads, self).run(folder)
