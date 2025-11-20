import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.abritamr.mainabritamr import main


class TestAbriTAMR(CamelTestSuite):
    """
    Tests the main AbriTAMR script.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_abritamr_standalone(self) -> None:
        """
        Tests the AbriTAMR standalone pipeline with fasta files.
        :return: None
        """
        path_report_html = self.running_dir / 'out' / 'report.html'
        path_report_tsv = self.running_dir / 'out' / 'report.tsv'

        result = cliutils.invoke(main, [
            '--fasta', str(self.input_fasta_file),
            '--output-html', str(path_report_html),
            '--output-dir', str(path_report_html.parent),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_report_tsv),
            '--threads', '2',
            '--species', 'Salmonella'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_html.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
