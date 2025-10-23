import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.utils import vcfutils
from camel.scripts.variantcalling.samtools.maincallingsamtools import MainCalling


class TestVariantCalling(CamelTestSuite):
    """
    Tests the variant calling tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('variant_calling')
    input_bam_file = test_file_dir / 'sars_cov_2-illumina.bam'
    input_bam_file_ont = test_file_dir / 'sars_cov_2-ont.bam'
    input_fasta_ref_file = test_file_dir / 'sars_cov_2-wuhan.fasta'

    def test_variant_calling(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'detected_variants.vcf'
        args = [
            '--bam', str(TestVariantCalling.input_bam_file),
            '--reference', str(TestVariantCalling.input_fasta_ref_file),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--ploidy', '1'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)

    def test_variant_calling_consensus(self) -> None:
        """
        Tests the variant calling main script to generate the consensus sequences.
        :return: None
        """
        output_file_vcf = self.running_dir / 'detected_variants.vcf'
        output_file_fasta = self.running_dir / 'consensus.fasta'
        args = [
            '--bam', str(TestVariantCalling.input_bam_file),
            '--reference', str(TestVariantCalling.input_fasta_ref_file),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--output-consensus', str(output_file_fasta),
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_fasta.stat().st_size, 0)

    def test_variant_calling_ont(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'detected_variants.vcf'
        args = [
            '--bam', str(TestVariantCalling.input_bam_file_ont),
            '--reference', str(TestVariantCalling.input_fasta_ref_file),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--ploidy', '1',
            '--input-type', 'ont'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)

if __name__ == '__main__':
    unittest.main()
