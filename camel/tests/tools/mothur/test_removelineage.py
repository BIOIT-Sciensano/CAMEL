import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurremovelineage import MothurRemoveLineage


class TestRemoveLineage(CamelTestSuite):
    """
    Tests Mothur remove.lineage.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'SAMPLE.trimmed_reads_1P.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'SAMPLE.trimmed_reads_1P.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table')
    FILE_TSV_TAXONOMY = ToolIOFile(test_file_dir / 'SAMPLE.trimmed_reads_1P.trim.contigs.good.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy')

    def test_removelineage(self) -> None:
        """
        Tests Mothur remove.lineage.
        :return: None
        """
        removelineage = MothurRemoveLineage()
        removelineage.add_input_files({
            'FASTA': [TestRemoveLineage.FILE_FASTA],
            'TSV_Counts': [TestRemoveLineage.FILE_TSV_COUNTS],
            'TSV_Taxonomy': [TestRemoveLineage.FILE_TSV_TAXONOMY]
        })
        removelineage.update_parameters(taxon='Chloroplast-Mitochondria-unknown-Archaea-Eukaryota')
        removelineage.run(self.running_dir)
        self.assertTrue('TSV_Taxonomy' in removelineage.tool_outputs, "No taxonomy output generated")
        self.assertTrue('TSV_Accnos' in removelineage.tool_outputs, "No accnos file generated")
        self.assertTrue('TSV_Counts' in removelineage.tool_outputs, "No counts file generated")
        self.assertTrue('FASTA' in removelineage.tool_outputs, "No fasta file generated")
        chimera_output = Path(removelineage.tool_outputs['TSV_Taxonomy'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(removelineage.tool_outputs['TSV_Accnos'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)
        chimera_output = Path(removelineage.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)
        accnos_output = Path(removelineage.tool_outputs['FASTA'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
