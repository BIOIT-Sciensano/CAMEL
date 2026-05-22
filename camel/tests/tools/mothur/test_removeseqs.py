import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurremoveseqs import MothurRemoveSeqs


class TestRemoveSeqs(CamelTestSuite):
    """
    Tests Mothur remove.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.fasta')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.count_table')
    FILE_TSV_ACCNOS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos')

    def test_removeseqs(self) -> None:
        """
        Tests Mothur remove.seqs.
        :return: None
        """
        removeseqs = MothurRemoveSeqs()
        removeseqs.add_input_files({
            'FASTA': [TestRemoveSeqs.FILE_FASTA],
            'TSV_Counts': [TestRemoveSeqs.FILE_TSV_COUNTS],
            'TSV_Accnos': [TestRemoveSeqs.FILE_TSV_ACCNOS]
        })
        removeseqs.run(self.running_dir)
        self.assertTrue('FASTA' in removeseqs.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Counts' in removeseqs.tool_outputs, "No counts file generated")
        chimera_output = Path(removeseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(removeseqs.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
