import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.utils import vcfutils
from camel.scripts.freebayes.mainfreebayes import MainFreebayesCalling


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
        args = [
            '--bam', str(self.FILE_BAM_ILLUMINA),
            '--reference', str(self.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output', str(output_file_vcf),
            '--ploidy', '1',
            '--standard-filters'
        ]
        main = MainFreebayesCalling(args)
        main.run()
        self.assertGreater(vcfutils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
