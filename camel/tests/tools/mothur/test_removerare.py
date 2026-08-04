import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurremoverare import MothurRemoveRare


class TestRemoveRare(CamelTestSuite):
    """
    Tests Mothur remove.rare.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')

    def test_removerare(self) -> None:
        """
        Tests Mothur remove.rare.
        :return: None
        """
        removerare = MothurRemoveRare()
        removerare.add_input_files({
            'TSV_List': [TestRemoveRare.FILE_TSV_LIST],
            'TSV_Counts': [TestRemoveRare.FILE_TSV_COUNTS]
        })
        removerare.update_parameters(nseqs=20)
        removerare.run(self.running_dir)
        self.assertTrue('TSV_List' in removerare.tool_outputs, "No list output generated")
        self.assertTrue('TSV_Counts' in removerare.tool_outputs, "No counts file generated")
        chimera_output = Path(removerare.tool_outputs['TSV_List'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(removerare.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
