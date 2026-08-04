import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import ToolExecutionError
from camel.app.tools.mothur.mothursummaryseqs import MothurSummarySeqs


class TestSummarySeqs(CamelTestSuite):
    """
    Tests Mothur summary.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.fasta')
    FILE_COUNT_TABLE = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.fake_group.count_table')

    def test_summaryseqs(self) -> None:
        """
        Tests Mothur summary.seqs.
        :return: None
        """
        summaryseqs = MothurSummarySeqs()
        summaryseqs.add_input_files({
            'FASTA': [TestSummarySeqs.FILE_FASTA]
        })
        summaryseqs.run(self.running_dir)
        self.assertTrue('TSV_Summary' in summaryseqs.tool_outputs, "No summapy output generated")
        self.assertTrue('TSV_Stats' in summaryseqs.tool_outputs, "No stats output generated")
        summary_output = Path(summaryseqs.tool_outputs['TSV_Summary'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        stats_output = Path(summaryseqs.tool_outputs['TSV_Stats'][0].path)
        self.assertTrue(stats_output.exists())
        self.assertGreater(stats_output.stat().st_size, 0)

    def test_summaryseqs_err_out(self) -> None:
        """
        Tests whether Mothur summary.seqs correctly raises a ToolExecutionError.
        :return: None
        """
        summaryseqs = MothurSummarySeqs()
        summaryseqs.add_input_files({
            'FASTA': [TestSummarySeqs.FILE_FASTA],
            'TSV_Counts': [TestSummarySeqs.FILE_COUNT_TABLE]
        })
        with self.assertRaises(ToolExecutionError):
            summaryseqs.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
