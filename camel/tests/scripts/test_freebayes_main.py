import unittest

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.freebayes.mainfreebayes import main


class TestFreebayesMain(CamelTestSuite):
    """
    Tests for the freebayes tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('clair3')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'bsubtilis_illumina.bam')

    def test_freebayes_maincalling(self) -> None:
        """
        Testing the freebayes main script.
        :return: None
        """
        output_file_vcf = self.running_dir / 'variants.vcf'
        result = cliutils.invoke(main, [
            '--bam', str(self.FILE_BAM_ILLUMINA),
            '--reference', str(self.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--ploidy', '1',
            '--standard-filters'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
