import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.checkv.maincheckv import main
from camel.tests import resourceIntensiveTest, longRunningTest


class TestCheckVMain(CamelTestSuite):
    """
    Tests the CheckV main script.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkv')
    input_fasta = test_file_dir / 'contigs_hev.fasta'

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkv_main_script(self) -> None:
        """
        Tests the CheckV main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        result = cliutils.invoke(main, [
            '--fasta', str(TestCheckVMain.input_fasta),
            '--working-dir', str(self.running_dir),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
