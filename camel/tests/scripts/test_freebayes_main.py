import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
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
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)


if __name__ == '__main__':
    unittest.main()
