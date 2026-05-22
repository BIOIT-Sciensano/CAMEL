import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurchimeravsearch import MothurChimeraVsearch


class TestChimeraVsearch(CamelTestSuite):
    """
    Tests Mothur chimera.vsearch.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.fasta')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.count_table')

    def test_chimeravsearch(self) -> None:
        """
        Tests Mothur chimera.vsearch.
        :return: None
        """
        chimeravsearch = MothurChimeraVsearch()
        chimeravsearch.add_input_files({
            'FASTA': [TestChimeraVsearch.FILE_FASTA],
            'TSV_Counts': [TestChimeraVsearch.FILE_TSV_COUNTS]
        })
        chimeravsearch.run(self.running_dir)
        self.assertTrue('TSV_Chimeras' in chimeravsearch.tool_outputs, "No chimeras output generated")
        self.assertTrue('TSV_Accnos' in chimeravsearch.tool_outputs, "No accnos file generated")
        self.assertTrue('FASTA' in chimeravsearch.tool_outputs, "No fasta output generated")
        self.assertTrue('TSV_Counts' in chimeravsearch.tool_outputs, "No counts file generated")
        chimera_output = Path(chimeravsearch.tool_outputs['TSV_Chimeras'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(chimeravsearch.tool_outputs['TSV_Accnos'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)
        fasta_output = Path(chimeravsearch.tool_outputs['FASTA'][0].path)
        self.assertTrue(fasta_output.exists())
        self.assertGreater(fasta_output.stat().st_size, 0)
        counts_output = Path(chimeravsearch.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(counts_output.exists())
        self.assertGreater(counts_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
