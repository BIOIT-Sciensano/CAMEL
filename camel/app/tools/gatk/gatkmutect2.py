import os

from camel.app.tools.gatk.gatk import GATK
from camel.app.core.io.tooliofile import ToolIOFile


class GATKMuTect2(GATK):
    """
    ==================
    GATK MuTect2 v 3.7
    ==================
    Performs variant calling.
    "MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller."
    web: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php

    Required inputs:
    ----------------
    "BAM_TUMOR":        ToolIOFile object. BAM file with tumour data

    Optional input:
    ---------------
    "BAM_NORMAL":       ToolIOFile object. BAM file with normal data for tumor-normal matching.
    "FASTA_REF":        ToolIOFile object. FASTA file containing the reference genome.
    "VCF_DBSNP":        ToolIOFile object. DbSNP reference vcf file location.
    "TXT_intervals":    ToolIOFile object. Intervals list to restrict search by GATK. Accelerates analysis. Bed or GATK intervals list 

    Output:
    -------
    "VCF":              ToolIOFile object. VCF file.

    Mandatory parameters:
    ---------------------
    - output_callstats_file     default value:  call_stats.txt
    - output_vcf_file           default value:  variant_call.vcf
    """

    def __init__(self):
        """
        Initialize GATKMuTect2 tool.
        :return: None
        """
        super(GATKMuTect2, self).__init__('gatk MuTect2', '3.7')

        self._required_inputs = ['BAM_TUMOR', 'FASTA_REF']
        self._specific_parameters.extend(["output_bam", "output_active_region_igv"])

    def _set_input(self):
        """
        Add the input specification in the input_string: 
        - reference fasta (in parent class)
        - intervals to work on (in parent class, optional)
        - DBSNP VCF file
        - tumour BAM file
        - normal tissue BAM file (optional)
        :return: None
        """
        super(GATKMuTect2, self)._set_input()

        # required tumor bam file
        self._input_string += "-I:tumor {} ".format(self._tool_inputs['BAM_TUMOR'][0].path)

        # optional normal bam file
        if 'BAM_NORMAL' in self._tool_inputs:
            self._input_string += "-I:normal {} ".format(self._tool_inputs['BAM_NORMAL'][0].path)

        # optional reference dbSNP vcf
        if 'VCF_DBSNP' in self._tool_inputs:
            self._input_string += "--dbsnp {} ".format(self._tool_inputs['VCF_DBSNP'][0].path)

        # optional reference cosmic vcf
        if 'VCF_COSMIC' in self._tool_inputs:
            self._input_string += "--cosmic {} ".format(self._tool_inputs['VCF_COSMIC'][0].path)

        # optional panel of normals vcf
        if 'VCF_PON' in self._tool_inputs:
            self._input_string += "--normal_panel {} ".format(self._tool_inputs['VCF_PON'][0].path)

    def _set_output(self):
        """
        Supersedes parent class function.
        Set the output specifications in the Camel output list: 
        - vcf file
        :return: None
        """
        self._tool_outputs['VCF'] = [
            ToolIOFile(self._folder / self.get_param_value('output_vcf_file'))]
        if 'output_bam' in self._parameters:
            if 'output_bam_file' not in self._params:
                self.update_parameters(output_bam_file=self.get_param('output_bam_file').value)
            self._tool_outputs['BAM'] = [ToolIOFile(self._folder / self.get_param_value('output_bam_file'))]
            self._tool_outputs['BAI'] = [
                ToolIOFile(self._folder / os.path.splitext(self.get_param_value('output_bam_file'))[0] + ".bai")]

    def _build_command(self):
        """
        Build the command to run tool. Supersedes that of parent class.
        :return: None
        """
        if "output_active_region_igv" in self._parameters:
            if "active_region_igv_file" not in self._parameters:
                self.update_parameters(active_region_igv_file=self.get_param("active_region_igv_file").value)

        self._option_string += " ".join(self._build_options(excluded_parameters=self._specific_parameters))

        self._command.command = " ".join([
            self._tool_command, self._input_string, self._output_string, self._option_string
        ])
