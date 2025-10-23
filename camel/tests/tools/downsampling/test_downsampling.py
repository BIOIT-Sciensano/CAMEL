import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.downsampling.downsamplecalculation import DownsampleCalculation
from camel.app.tools.pipelines.downsampling.fastqstats import FastqStats


class TestDownsampling(CamelTestSuite):
    """
    Tests the downsampling functionality.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('downsampling')
    FILE_FQ_FWD = ToolIOFile(test_file_dir / 'ecoli_1.fastq')
    FILE_FQ_REV = ToolIOFile(test_file_dir / 'ecoli_2.fastq')
    FILE_FQ_GZ_FWD = ToolIOFile(test_file_dir / 'ecoli_1.fastq.gz')
    FILE_FQ_GZ_REV = ToolIOFile(test_file_dir / 'ecoli_2.fastq.gz')

    def test_fastq_stats(self) -> None:
        """
        Tests the FASTQ stats tool.
        :return: None
        """
        fqtools_stats = FastqStats()
        fqtools_stats.add_input_files({
            'FASTQ': [TestDownsampling.FILE_FQ_FWD, TestDownsampling.FILE_FQ_REV]})
        fqtools_stats.run(self.running_dir)
        self.assertIn('stats', fqtools_stats.informs)
        self.assertEqual(len(fqtools_stats.informs['stats']), 2)
        self.assertIn('nb_of_sequences', fqtools_stats.informs['stats'][0])
        self.assertIn('nb_of_bases', fqtools_stats.informs['stats'][0])

    def test_fastq_stats_gzip(self) -> None:
        """
        Tests the FASTQ stats tool with gzipped input.
        :return: None
        """
        fqtools_stats = FastqStats()
        fqtools_stats.add_input_files({
            'FASTQ': [TestDownsampling.FILE_FQ_GZ_FWD, TestDownsampling.FILE_FQ_GZ_REV]})
        fqtools_stats.run(self.running_dir)
        self.assertIn('stats', fqtools_stats.informs)
        self.assertEqual(len(fqtools_stats.informs['stats']), 2)
        self.assertIn('nb_of_sequences', fqtools_stats.informs['stats'][0])
        self.assertIn('nb_of_bases', fqtools_stats.informs['stats'][0])

    def test_downsampling_calculate(self) -> None:
        """
        Tests the calculation of the downsampling statistics.
        """
        # Calculate FASTQ stats
        fqtools_stats = FastqStats()
        fqtools_stats.add_input_files({
            'FASTQ': [TestDownsampling.FILE_FQ_GZ_FWD, TestDownsampling.FILE_FQ_GZ_REV]})
        fqtools_stats.run(self.running_dir)

        # Calculate downsampling stats
        ds_calc = DownsampleCalculation()
        ds_calc.add_input_informs({'stats': fqtools_stats.informs})
        ds_calc.update_parameters(size_ref_genome=10_000, cov_target=1, is_paired=True)
        ds_calc.run(self.running_dir)
        self.verify_output_files(ds_calc, 'JSON')
        self.assertIn('stats', ds_calc.informs)

    def test_downsampling_calculate_se(self) -> None:
        """
        Tests the calculation of the downsampling statistics.
        """
        # Calculate FASTQ stats
        fqtools_stats = FastqStats()
        fqtools_stats.add_input_files({'FASTQ': [TestDownsampling.FILE_FQ_GZ_FWD]})
        fqtools_stats.run(self.running_dir)

        # Calculate downsampling stats
        ds_calc = DownsampleCalculation()
        ds_calc.add_input_informs({'stats': fqtools_stats.informs})
        ds_calc.update_parameters(size_ref_genome=10_000, cov_target=1, is_paired=False)
        ds_calc.run(self.running_dir)
        self.verify_output_files(ds_calc, 'JSON')
        self.assertIn('stats', ds_calc.informs)


if __name__ == '__main__':
    unittest.main()
