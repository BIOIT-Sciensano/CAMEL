import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.spatyping.mainspatyping import MainSpaTyping


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
        args = [
            '--fasta', str(TestSpaTyping.input_fasta_file),
            '--input-type', 'fasta',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-path', '/db/pipelines/saureus/spatyping'
        ]
        main = MainSpaTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
