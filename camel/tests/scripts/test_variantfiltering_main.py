import unittest

from camelcore.app.utils import vcfutils

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.variantcalling.mainfilteringsamtools import main


class TestVariantFilteringMain(CamelTestSuite):
    """
    Tests the variant filters main script.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('variant_calling')

    PATH_VCF_GZ_UNFILTERED = test_file_dir / 'unfiltered_variants-myco.vcf.gz'
    PATH_BAM = test_file_dir / 'alignment.bam'
    PATH_BED = test_file_dir / 'h37Rv.bed'

    def test_variant_filtering_main(self) -> None:
        """
        Tests the main script for the variant filtering.
        :return: None
        """
        number_variants_in = vcfutils.count_variants(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED)
        output_file_vcf = self.running_dir / 'filtered_variants.vcf'
        output_file_stats = self.running_dir / 'filter_stats.txt'
        result = cliutils.invoke(main, [
            '--vcf', str(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED),
            '--bam', str(TestVariantFilteringMain.PATH_BAM),
            '--working-dir', str(self.running_dir),
            '--output-vcf', str(output_file_vcf),
            '--output-stats', str(output_file_stats)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(output_file_stats.stat().st_size, 0)
        self.assertLess(vcfutils.count_variants(output_file_vcf), number_variants_in)

    def test_variant_filtering_main_soft(self) -> None:
        """
        Tests the main script for the variant filtering with soft masking.
        :return: None
        """
        number_variants_in = vcfutils.count_variants(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED)
        output_file_vcf = self.running_dir / 'filtered_variants.vcf'
        output_file_stats = self.running_dir / 'filter_stats.txt'
        result = cliutils.invoke(main, [
            '--vcf', str(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED),
            '--bam', str(TestVariantFilteringMain.PATH_BAM),
            '--working-dir', str(self.running_dir),
            '--output-vcf', str(output_file_vcf),
            '--output-stats', str(output_file_stats),
            '--soft-filter'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(output_file_stats.stat().st_size, 0)
        self.assertLess(vcfutils.count_variants(output_file_vcf), number_variants_in)

    def test_variant_filtering_main_no_bam(self) -> None:
        """
        Tests the main script for the variant filtering.
        :return: None
        """
        number_variants_in = vcfutils.count_variants(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED)
        output_file_vcf = self.running_dir / 'filtered_variants.vcf'
        result = cliutils.invoke(main, [
            '--vcf', str(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED),
            '--bam', str(TestVariantFilteringMain.PATH_BAM),
            '--working-dir', str(self.running_dir),
            '--output-vcf', str(output_file_vcf),
            '--min-snp-quality', '30'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertLess(vcfutils.count_variants(output_file_vcf), number_variants_in)

    def test_variant_filtering_main_bed_input(self) -> None:
        """
        Tests the main script for the variant filtering with a BED file to filter problematic regions.
        :return: None
        """
        number_variants_in = vcfutils.count_variants(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED)
        output_file_vcf = self.running_dir / 'filtered_variants.vcf'
        result = cliutils.invoke(main, [
            '--vcf', str(TestVariantFilteringMain.PATH_VCF_GZ_UNFILTERED),
            '--bam', str(TestVariantFilteringMain.PATH_BAM),
            '--bed', str(TestVariantFilteringMain.PATH_BED),
            '--working-dir', str(self.running_dir),
            '--output-vcf', str(output_file_vcf),
            '--min-snp-quality', '30'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertLess(vcfutils.count_variants(output_file_vcf), number_variants_in)


if __name__ == '__main__':
    unittest.main()
