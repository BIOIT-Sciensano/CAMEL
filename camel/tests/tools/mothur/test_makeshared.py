import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurmakeshared import MothurMakeShared


class TestMakeShared(CamelTestSuite):
    """
    Tests Mothur make.shared.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')

    def test_makeshared(self) -> None:
        """
        Tests Mothur make.shared.
        :return: None
        """
        makeshared = MothurMakeShared()
        makeshared.add_input_files({
            'TSV_List': [TestMakeShared.FILE_TSV_LIST],
            'TSV_Counts': [TestMakeShared.FILE_TSV_COUNTS]
        })
        makeshared.run(self.running_dir)
        self.assertTrue('TSV_Shared' in makeshared.tool_outputs, "No shared output generated")
        chimera_output = Path(makeshared.tool_outputs['TSV_Shared'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
