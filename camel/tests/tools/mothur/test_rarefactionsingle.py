import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurrarefactionsingle import MothurRarefactionSingle


class TestRarefactionSingle(CamelTestSuite):
    """
    Tests Mothur rarefaction.single.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list')

    def test_rarefactionsingle(self) -> None:
        """
        Tests Mothur rarefaction.single.
        :return: None
        """
        rarefactionsingle = MothurRarefactionSingle()
        rarefactionsingle.add_input_files({
            'TSV_List': [TestRarefactionSingle.FILE_TSV_LIST]
        })
        rarefactionsingle.run(self.running_dir)
        self.assertTrue('TSV' in rarefactionsingle.tool_outputs, "No rarefaction output generated")
        chimera_output = Path(rarefactionsingle.tool_outputs['TSV'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
