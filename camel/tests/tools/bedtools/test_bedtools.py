import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.bedtools.bedtoolsgenomecov import BedtoolsGenomecov


class TestBedtools(CamelTestSuite):
    """
    Tests the bedtools tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('bedtools')
    TEST_FILE_BAM = ToolIOFile(test_file_dir / 'samtools_sort.bam')

    def test_bedtools_genomecov_bed_out(self) -> None:
        """
        Tests bedtools genomecov.
        :return: None
        """
        bedtools_genome_cov = BedtoolsGenomecov()
        bedtools_genome_cov.add_input_files({'BAM': [TestBedtools.TEST_FILE_BAM]})
        bedtools_genome_cov.update_parameters(Depth=False, BedGraphWithZeroCoverage=True)
        bedtools_genome_cov.run(self.running_dir)
        self.verify_output_files(bedtools_genome_cov, 'BED')

    def test_bedtools_genomecov_tsv_out(self) -> None:
        """
        Tests bedtools genomecov with tabular output format.
        :return: None
        """
        bedtools_genome_cov = BedtoolsGenomecov()
        bedtools_genome_cov.add_input_files({'BAM': [TestBedtools.TEST_FILE_BAM]})
        bedtools_genome_cov.update_parameters(Depth=True)
        bedtools_genome_cov.run(self.running_dir)
        self.verify_output_files(bedtools_genome_cov, 'TSV')

    def test_bedtools_genome_cov_invalid_input(self) -> None:
        """
        Tests if bedtools genomecov with invalid input raises an exception.
        :return: None
        """
        with self.assertRaises(InvalidToolInputError):
            bedtools_genome_cov = BedtoolsGenomecov()
            bedtools_genome_cov.add_input_files({'TXT': [TestBedtools.TEST_FILE_BAM]})
            bedtools_genome_cov.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
