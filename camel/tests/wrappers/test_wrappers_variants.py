import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.scriptutils.fastqinput import FastqInput
from camel.app.wrappers.variantcallingwrapper import VariantCallingWrapper
from camel.app.wrappers.variantfilteringwrapper import VariantFilteringWrapper


class TestWrappersVariants(CamelTestSuite):
    """
    Tests the Snakemake variant calling related wrappers.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir()

    # Variant calling
    input_fastq_pe = [
        test_file_dir / 'variant_calling' / 'toy' / 'r_1.fq',
        test_file_dir / 'variant_calling' / 'toy' / 'r_2.fq']
    input_fasta_ref = test_file_dir / 'variant_calling' / 'toy' / 'toy.fasta'

    # Variant filtering
    input_vcf = test_file_dir / 'variant_calling' / 'unfiltered_variants-myco.vcf.gz'
    input_bam = test_file_dir / 'variant_calling' / 'alignment.bam'

    def test_variant_calling_illumina(self) -> None:
        """
        Tests the variant calling wrapper.
        :return: None
        """
        wrapper = VariantCallingWrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWrappersVariants.input_fastq_pe])
        ref_info = {
            'name': TestWrappersVariants.input_fasta_ref.name,
            'path': str(TestWrappersVariants.input_fasta_ref)
        }
        wrapper.run(ref_info, 'test_sample', fastq_input, 'illumina',{}, 1)
        self.assertGreater(wrapper.output.bam_file.size, 0)
        self.assertGreater(wrapper.output.vcf_unfiltered.size, 0)

    def test_variant_filtering(self) -> None:
        """
        Tests the variant filtering wrapper.
        :return: None
        """
        wrapper = VariantFilteringWrapper(self.running_dir)
        wrapper.run(
            'test_sample', TestWrappersVariants.input_vcf, TestWrappersVariants.input_bam, 'illumina', {}, 1)
        self.assertGreater(Path(wrapper.output.vcf_filtered.path).stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
