import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.variantcalling.clair3.maincallingclair3 import MainCalling


class TestClair3(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('clair3')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'bsubtilis_illumina.bam')
    FILE_BAM_ONT = ToolIOFile(test_file_dir / 'bsubtilis_ont.bam')

    def test_variant_calling_ont(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'clair3_output.vcf'
        args = [
            '--bam', str(TestClair3.FILE_BAM_ONT),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', '/db/clair3/models/ont/',
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ont'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)

    def test_variant_calling_illumina(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'clair3_output.vcf'
        args = [
            '--bam', str(TestClair3.FILE_BAM_ILLUMINA),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', '/db/clair3/models/ilmn/',
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ilmn'
        ]
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
