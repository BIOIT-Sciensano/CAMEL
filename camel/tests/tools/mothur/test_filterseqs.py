import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurfilterseqs import MothurFilterSeqs


class TestFilterSeqs(CamelTestSuite):
    """
    Tests Mothur filter.seqs.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.align')

    def test_filterseqs(self) -> None:
        """
        Tests Mothur filter.seqs.
        :return: None
        """
        filterseqs = MothurFilterSeqs()
        filterseqs.add_input_files({
            'FASTA': [TestFilterSeqs.FILE_FASTA]
        })
        filterseqs.run(self.running_dir)
        self.assertTrue('FASTA' in filterseqs.tool_outputs, "No FASTA output generated")
        summary_output = Path(filterseqs.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
