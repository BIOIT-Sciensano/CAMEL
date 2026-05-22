import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurgetoturep import MothurGetOturep


class TestGetOturep(CamelTestSuite):
    """
    Tests Mothur get.oturep.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta')

    def test_getoturep(self) -> None:
        """
        Tests Mothur get.oturep.
        :return: None
        """
        getoturep = MothurGetOturep()
        getoturep.add_input_files({
            'TSV_List': [TestGetOturep.FILE_TSV_LIST],
            'TSV_Counts': [TestGetOturep.FILE_TSV_COUNTS],
            'FASTA': [TestGetOturep.FILE_FASTA]
        })
        getoturep.update_parameters(method='abundance')
        getoturep.run(self.running_dir)
        self.assertTrue('FASTA' in getoturep.tool_outputs, "No fasta output generated")
        self.assertTrue('TSV_Counts' in getoturep.tool_outputs, "No counts file generated")
        chimera_output = Path(getoturep.tool_outputs['FASTA'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(getoturep.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
