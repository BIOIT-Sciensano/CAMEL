import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurcluster import MothurCluster


class TestCluster(CamelTestSuite):
    """
    Tests Mothur cluster.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_DIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')

    def test_cluster(self) -> None:
        """
        Tests Mothur cluster.
        :return: None
        """
        cluster = MothurCluster()
        cluster.add_input_files({
            'DIST': [TestCluster.FILE_DIST],
            'TSV_Counts': [TestCluster.FILE_TSV_COUNTS]
        })
        cluster.update_parameters(method='opti')
        cluster.run(self.running_dir)
        self.assertTrue('TSV_List' in cluster.tool_outputs, "No list output generated")
        chimera_output = Path(cluster.tool_outputs['TSV_List'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
