import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothuruniqueseqs import MothurUniqueSeqs


class TestUniqueSeqs(CamelTestSuite):
    """
    Tests Mothur unique.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.fasta')

    def test_uniqueseqs(self) -> None:
        """
        Tests Mothur unique.seqs.
        :return: None
        """
        uniqueseqs = MothurUniqueSeqs()
        uniqueseqs.add_input_files({
            'FASTA': [TestUniqueSeqs.FILE_FASTA]
        })
        uniqueseqs.run(self.running_dir)
        self.assertTrue('FASTA' in uniqueseqs.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Counts' in uniqueseqs.tool_outputs, "No file with counts generated")
        summary_output = Path(uniqueseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        stats_output = Path(uniqueseqs.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(stats_output.exists())
        self.assertGreater(stats_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
