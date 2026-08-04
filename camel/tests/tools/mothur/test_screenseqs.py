import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurscreenseqs import MothurScreenSeqs


class TestScreenSeqs(CamelTestSuite):
    """
    Tests Mothur screen.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.fasta')

    def test_screenseqs(self) -> None:
        """
        Tests Mothur screen.seqs.
        :return: None
        """
        screenseqs = MothurScreenSeqs()
        screenseqs.add_input_files({
            'FASTA': [TestScreenSeqs.FILE_FASTA]
        })
        screenseqs.update_parameters(maxambig=15, maxhomop=5)
        screenseqs.run(self.running_dir)
        self.assertTrue('FASTA' in screenseqs.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Bad' in screenseqs.tool_outputs, "No file with bad sequence IDs generated")
        summary_output = Path(screenseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        stats_output = Path(screenseqs.tool_outputs['TSV_Bad'][0].path)
        self.assertTrue(stats_output.exists())
        self.assertGreater(stats_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
