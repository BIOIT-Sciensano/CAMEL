import unittest
from pathlib import Path

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.utils import vcfutils
from camel.scripts.variantcalling.maincallingclair3 import main


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
        result = cliutils.invoke(main, [
            '--bam', str(TestClair3.FILE_BAM_ONT),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ont'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)

    def test_variant_calling_illumina(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'clair3_output.vcf'
        result = cliutils.invoke(main, [
            '--bam', str(TestClair3.FILE_BAM_ILLUMINA),
            '--reference', str(TestClair3.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--model-path', str(Path(config.dir_db, 'clair3', 'models', 'ilmn')),
            '--haploid-precise',
            '--no-phasing',
            '--include-ctgs',
            '--platform', 'ilmn'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_vcf.stat().st_size, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
