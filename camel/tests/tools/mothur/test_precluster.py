import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurprecluster import MothurPreCluster


class TestPreCluster(CamelTestSuite):
    """
    Tests Mothur pre.cluster.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.filter.unique.fasta')
    FILE_COUNT_TABLE = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.fake_group.count_table')
    FILE_COUNT_TABLE_NO_GROUPS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.count_table')

    def test_precluster(self) -> None:
        """
        Tests Mothur pre.cluster.
        :return: None
        """
        precluster = MothurPreCluster()
        precluster.add_input_files({
            'FASTA': [TestPreCluster.FILE_FASTA],
            'TSV_Counts': [TestPreCluster.FILE_COUNT_TABLE]
        })
        precluster.update_parameters(diffs=2, processors=6)
        precluster.run(self.running_dir)
        self.assertTrue('FASTA' in precluster.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Counts' in precluster.tool_outputs, "No counts output generated")
        summary_output = Path(precluster.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        summary_output = Path(precluster.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)

    def test_precluster_counts_without_groups(self) -> None:
        """
        Tests Mothur pre.cluster.
        :return: None
        """
        precluster = MothurPreCluster()
        precluster.add_input_files({
            'FASTA': [TestPreCluster.FILE_FASTA],
            'TSV_Counts': [TestPreCluster.FILE_COUNT_TABLE_NO_GROUPS]
        })
        precluster.update_parameters(diffs=2)
        precluster.run(self.running_dir)
        self.assertTrue('FASTA' in precluster.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Counts' in precluster.tool_outputs, "No counts output generated")
        summary_output = Path(precluster.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        summary_output = Path(precluster.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
