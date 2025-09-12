import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.checkm.maincheckm import MainCheckM
from camel.tests import longRunningTest, resourceIntensiveTest


class TestCheckMMain(CamelTestSuite):
    """
    Tests the CheckM main script.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkm')
    input_fasta = test_file_dir / 'contigs_neisseria.fasta'

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkm_main_script(self) -> None:
        """
        Tests the CheckM main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        checkm_main = MainCheckM([
            '--fasta', str(TestCheckMMain.input_fasta), TestCheckMMain.input_fasta.name,
            '--working-dir', str(self.running_dir),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--reduced_tree'
        ])
        checkm_main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkm_main_script_spaces_in_name(self) -> None:
        """
        Tests the CheckM main script, with spaces in the input file name (for Galaxy).
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        checkm_main = MainCheckM([
            '--fasta', str(TestCheckMMain.input_fasta), '"SPAdes on data 126 and data 125 - Contigs"',
            '--working-dir', str(self.running_dir),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--reduced_tree'
        ])
        checkm_main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
