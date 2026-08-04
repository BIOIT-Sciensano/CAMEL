import unittest
from pathlib import Path

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.spatyping.mainspatyping import main


class TestSpaTyping(CamelTestSuite):
    """
    Tests the spa typing tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('spatyping')
    input_fasta_file = test_file_dir / 'saureus_ref.fasta'

    def test_spatyping_fasta(self) -> None:
        """
        Tests the spa typing main script
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        result = cliutils.invoke(main, [
            '--fasta', str(TestSpaTyping.input_fasta_file),
            '--input-type', 'fasta',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db', str(Path(config.dir_db, 'pipelines', 'saureus', 'spatyping'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
