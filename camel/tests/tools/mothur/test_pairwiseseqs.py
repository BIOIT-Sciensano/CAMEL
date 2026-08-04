import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurpairwiseseqs import MothurPairwiseSeqs


class TestPairwiseSeqs(CamelTestSuite):
    """
    Tests Mothur pairwise.seqs.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.subsample.fasta')

    def test_pairwiseseqs(self) -> None:
        """
        Tests Mothur pairwise.seqs.
        :return: None
        """
        pairwiseseqs = MothurPairwiseSeqs()
        pairwiseseqs.add_input_files({
            'FASTA': [TestPairwiseSeqs.FILE_FASTA]
        })
        pairwiseseqs.update_parameters(processors=6)
        pairwiseseqs.run(self.running_dir)
        self.assertTrue('DIST' in pairwiseseqs.tool_outputs, "No distance matrix output generated")
        summary_output = Path(pairwiseseqs.tool_outputs['DIST'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
