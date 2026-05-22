import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.amrfinder.mainamrfinder import main


class TestAMRFinder(CamelTestSuite):
    """
    Tests the AMRFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('amrfinder')
    FILE_FASTA = ToolIOFile(test_file_dir / 'test_dna.fa')
    DIR_DB = config.dir_db / 'amrfinder' / 'v4' / 'latest'

    def test_amrfinder_main(self) -> None:
        """
        Tests the AMRFinder main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        result = cliutils.invoke(
            main, [
                '--fasta', str(TestAMRFinder.FILE_FASTA),
                '--db', str(TestAMRFinder.DIR_DB),
                '--output-html', str(path_report_out),
                '--output-dir', str(path_report_out.parent),
                '--working-dir', str(self.running_dir),
            ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_amrfinder_main_params(self) -> None:
        """
        Tests the AMRFinder main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        result = cliutils.invoke(
            main, [
            '--fasta', str(TestAMRFinder.FILE_FASTA),
            '--fasta-name', 'assembly name with spaces.fasta',
            '--db', str(TestAMRFinder.DIR_DB),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--min-id', '95',
            '--min-cov', '60',
            '--organism', 'Escherichia'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_amrfinder_main_tsv_out(self) -> None:
        """
        Tests the AMRFinder main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        path_tsv_out = self.running_dir / 'extra_output_amrfinder.tsv'
        result = cliutils.invoke(
            main, [
            '--fasta', str(TestAMRFinder.FILE_FASTA),
            '--db', str(TestAMRFinder.DIR_DB),
            '--output-html', str(path_report_out),
            '--output-tsv', str(path_tsv_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_tsv_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
