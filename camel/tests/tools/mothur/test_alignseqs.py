import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothuralignseqs import MothurAlignSeqs


class TestAlignSeqs(CamelTestSuite):
    """
    Tests Mothur align.seqs.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'silva.v4.fasta')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.fasta')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.good.filter.unique.precluster.count_table')

    def test_alignseqs(self) -> None:
        """
        Tests Mothur align.seqs.
        :return: None
        """
        alignseqs = MothurAlignSeqs()
        alignseqs.add_input_files({
            'FASTA_Ref': [TestAlignSeqs.FILE_FASTA_REF],
            'FASTA': [TestAlignSeqs.FILE_FASTA]
        })
        alignseqs.update_parameters(processors=6)
        alignseqs.run(self.running_dir)
        self.assertTrue('FASTA' in alignseqs.tool_outputs, "No FASTA output generated")
        self.assertTrue('TSV_Report' in alignseqs.tool_outputs, "No report output generated")
        summary_output = Path(alignseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)
        summary_output = Path(alignseqs.tool_outputs['TSV_Report'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)

    def test_alignseqs_incomplete_input(self) -> None:
        """
        Tests that mothur_alignseqs fails with incomplete input.
        :return: None
        """
        alignseqs = MothurAlignSeqs()
        alignseqs.add_input_files({
            'FASTA_Ref': [TestAlignSeqs.FILE_FASTA_REF],
        })
        with self.assertRaises(InvalidToolInputError):
            alignseqs.run(self.running_dir)

    def test_alignseqs_incorrect_input(self) -> None:
        """
        Tests that mothur_alignseqs fails with incorrect (not allowed) input.
        :return: None
        """
        alignseqs = MothurAlignSeqs()
        alignseqs.add_input_files({
            'FASTA_Ref': [TestAlignSeqs.FILE_FASTA_REF],
            'TSV_Counts': [TestAlignSeqs.FILE_FASTA]
        })
        with self.assertRaises(InvalidToolInputError):
            alignseqs.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
