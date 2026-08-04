import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothursubsample import MothurSubSample


class TestSubSample(CamelTestSuite):
    """
    Tests Mothur sub.sample.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTA = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.trim.contigs.good.unique.fasta')

    def test_subsample(self) -> None:
        """
        Tests Mothur sub.sample.
        :return: None
        """
        subsample = MothurSubSample()
        subsample.add_input_files({
            'FASTA': [TestSubSample.FILE_FASTA]
        })
        subsample.run(self.running_dir)
        self.assertTrue('FASTA' in subsample.tool_outputs, "No FASTA output generated")
        summary_output = Path(subsample.tool_outputs['FASTA'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
