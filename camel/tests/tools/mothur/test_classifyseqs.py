import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurclassifyseqs import MothurClassifySeqs


class TestClassifySeqs(CamelTestSuite):
    """
    Tests Mothur classify.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'silva.pcr.align')
    FILE_TSV_TAXONOMY = ToolIOFile(test_file_dir / 'silva.tax')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')

    def test_classifyseqs(self) -> None:
        """
        Tests Mothur classify.seqs.
        :return: None
        """
        classifyseqs = MothurClassifySeqs()
        classifyseqs.add_input_files({
            'FASTA': [TestClassifySeqs.FILE_FASTA],
            'TSV_Counts': [TestClassifySeqs.FILE_TSV_COUNTS],
            'TSV_Taxonomy': [TestClassifySeqs.FILE_TSV_TAXONOMY],
            'FASTA_Ref': [TestClassifySeqs.FILE_FASTA_REF]
        })
        classifyseqs.update_parameters(cutoff=80, processors=8)
        classifyseqs.run(self.running_dir)
        self.assertTrue('TSV_Taxonomy' in classifyseqs.tool_outputs, "No taxonomy output generated")
        self.assertTrue('TSV_Summary' in classifyseqs.tool_outputs, "No tax summary file generated")
        chimera_output = Path(classifyseqs.tool_outputs['TSV_Taxonomy'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(classifyseqs.tool_outputs['TSV_Summary'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
