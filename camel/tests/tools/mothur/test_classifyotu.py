import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurclassifyotu import MothurClassifyOtu


class TestClassifyOtu(CamelTestSuite):
    """
    Tests Mothur classify.otu.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list')
    FILE_TSV_TAXONOMY = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')

    def test_classifyotu(self) -> None:
        """
        Tests Mothur classify.otu.
        :return: None
        """
        classifyotu = MothurClassifyOtu()
        classifyotu.add_input_files({
            'TSV_List': [TestClassifyOtu.FILE_TSV_LIST],
            'TSV_Counts': [TestClassifyOtu.FILE_TSV_COUNTS],
            'TSV_Taxonomy': [TestClassifyOtu.FILE_TSV_TAXONOMY]
        })
        classifyotu.run(self.running_dir)
        self.assertTrue('TSV_Taxonomy' in classifyotu.tool_outputs, "No taxonomy output generated")
        self.assertTrue('TSV_Summary' in classifyotu.tool_outputs, "No summary file generated")
        chimera_output = Path(classifyotu.tool_outputs['TSV_Taxonomy'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(classifyotu.tool_outputs['TSV_Summary'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
