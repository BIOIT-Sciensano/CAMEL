import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurcountseqs import MothurCountSeqs


class TestAlignSeqs(CamelTestSuite):
    """
    Tests Mothur align.seqs.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_NAMES = ToolIOFile(test_file_dir / 'example.name')

    def test_countseqs(self) -> None:
        """
        Tests Mothur align.seqs.
        :return: None
        """
        countseqs = MothurCountSeqs()
        countseqs.add_input_files({
            'TSV_Names': [TestAlignSeqs.FILE_TSV_NAMES]
        })
        countseqs.run(self.running_dir)
        self.assertTrue('TSV_Counts' in countseqs.tool_outputs, "No TSV counts output generated")
        summary_output = Path(countseqs.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(summary_output.exists())
        self.assertGreater(summary_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
