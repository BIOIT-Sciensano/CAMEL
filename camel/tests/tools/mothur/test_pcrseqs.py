import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurpcrseqs import MothurPcrSeqs


class TestPcrSeqs(CamelTestSuite):
    """
    Tests Mothur pcr.seqs.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'silva.bacteria.fasta')

    def test_pcrseqs(self) -> None:
        """
        Tests Mothur pcr.seqs.
        :return: None
        """
        pcrseqs = MothurPcrSeqs()
        pcrseqs.add_input_files({
            'FASTA': [TestPcrSeqs.FILE_FASTA]
        })
        pcrseqs.update_parameters(start=11895, end=25318, keepdots='F', processors=6)
        pcrseqs.run(self.running_dir)
        self.assertTrue('FASTA' in pcrseqs.tool_outputs, "No FASTA output generated")
        summary_output = Path(pcrseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
