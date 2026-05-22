import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.integronfinder.mainintegronfinder import main


class TestIntegronFinder(CamelTestSuite):
    """
    Tests for the IntegronFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('integron_finder')
    FILE_FASTA = ToolIOFile(test_file_dir / 'GCF_002202175.1.fasta')
    FILE_FASTA_E_COLI = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')

    def test_integron_finder_main_script(self) -> None:
        """
        Tests the IntegronFinder main script.
        """
        path_out_html = self.running_dir / 'out' / 'report.html'
        path_out_tsv = self.running_dir / 'out' / 'integrons.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestIntegronFinder.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_out_tsv),
            '--output-html', str(path_out_html),
            '--output-dir', str(path_out_html.parent),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)

    def test_integron_finder_main_script_ecoli(self) -> None:
        """
        Tests the IntegronFinder main script on E. coli.
        """
        path_out_html = self.running_dir / 'out' / 'report.html'
        path_out_tsv = self.running_dir / 'out' / 'integrons.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestIntegronFinder.FILE_FASTA_E_COLI),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_out_tsv),
            '--output-html', str(path_out_html),
            '--output-dir', str(path_out_html.parent),
            '--threads', '4'
        ])
        self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()
