from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.tools.gatk.gatk import GATK


class GATKMuTect2(GATK):
    """
    Calls somatic SNPs and indels, combining the somatic genotyping engine of MuTect with the assembly-based machinery
    of HaplotypeCaller.

    Required inputs:
    - BAM_TUMOR: Alignment file of the tumor sample
    - FASTA_REF: Reference genome

    Optional inputs:
    - BAM_NORMAL: Alignment file of the matched normal sample
    - VCF_DBSNP: dbSNP variants
    - VCF_COSMIC: COSMIC variants
    - VCF_PON: Panel of normals
    - TXT_intervals: Intervals to restrict the analysis to

    Outputs:
    - VCF: Variants
    - BAM / BAI: Assembled haplotypes (only if 'output_bam' or 'output_bam_file' is set)
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk MuTect2')

    def _check_parameters(self) -> None:
        """
        Checks if the parameters are valid.
        :return: None
        """
        super()._check_parameters()
        # The flags enable the corresponding output file parameters using their default file name
        if self.get_param_value('output_bam') and self.get_param_value('output_bam_file') is None:
            self.update_parameters(output_bam_file=self.get_param('output_bam_file').value)
        if self.get_param_value('output_active_region_igv') and self.get_param_value('active_region_igv_file') is None:
            self.update_parameters(active_region_igv_file=self.get_param('active_region_igv_file').value)

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['BAM_TUMOR', 'FASTA_REF'])
        super()._check_input()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [*super()._get_input_options(), f"-I:tumor {self._tool_inputs['BAM_TUMOR'][0].path}"]
        for key, option in [
                ('BAM_NORMAL', '-I:normal'),
                ('VCF_DBSNP', '--dbsnp'),
                ('VCF_COSMIC', '--cosmic'),
                ('VCF_PON', '--normal_panel')]:
            if key in self._tool_inputs:
                options.append(f'{option} {self._tool_inputs[key][0].path}')
        return options

    def _get_parameter_options(self) -> list[str]:
        """
        Returns the command line options for the parameters.
        :return: Parameter options
        """
        return self._build_options(excluded_parameters=['output_bam', 'output_active_region_igv'])

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(self.folder / self.get_param_value('output_vcf_file'))]
        if self.get_param_value('output_bam_file') is not None:
            path_bam = self.folder / self.get_param_value('output_bam_file')
            self._tool_outputs['BAM'] = [ToolIOFile(path_bam)]
            self._tool_outputs['BAI'] = [ToolIOFile(path_bam.with_suffix('.bai'))]
