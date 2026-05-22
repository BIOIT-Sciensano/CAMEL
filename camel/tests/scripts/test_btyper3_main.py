import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.btyper.mainbtyper import main


class TestBTyper3Main(CamelTestSuite):
    """
    Tests for BTyper main script.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('btyper')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bacillus_contigs.fasta')

    def test_btyper_main(self) -> None:
        """
        Tests the BTyper main script.
        :return: None
        """
        path_out_html = self.running_dir / 'out' / 'report.html'
        path_out_html.parent.mkdir()
        path_out_tsv = self.running_dir / 'out' / 'tabular.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestBTyper3Main.FILE_FASTA),
            '--fasta-name', str(TestBTyper3Main.FILE_FASTA),
            '--output-html', str(path_out_html),
            '--output-dir', str(path_out_html.parent),
            '--output-tsv', str(path_out_tsv),
            '--working-dir', str(self.running_dir),
            '--mlst', '--panc', '--bt', '--virulence'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertTrue(path_out_html.exists())
        self.assertTrue(path_out_tsv.exists())


if __name__ == '__main__':
    unittest.main()
