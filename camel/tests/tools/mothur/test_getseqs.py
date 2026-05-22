import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurgetseqs import MothurGetSeqs


class TestGetSeqs(CamelTestSuite):
    """
    Tests Mothur get.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_ACCNOS = ToolIOFile(test_file_dir / 'ITS.accnos')
    FILE_TSV_COUNTS = ToolIOFile(test_file_dir / 'ITS.count_table')

    def test_getseqs(self) -> None:
        """
        Tests Mothur get.seqs.
        :return: None
        """
        getseqs = MothurGetSeqs()
        getseqs.add_input_files({
            'TSV_Accnos': [TestGetSeqs.FILE_TSV_ACCNOS],
            'TSV_Counts': [TestGetSeqs.FILE_TSV_COUNTS]
        })
        getseqs.run(self.running_dir)
        self.assertTrue('TSV_Counts' in getseqs.tool_outputs, "No counts output generated")
        chimera_output = Path(getseqs.tool_outputs['TSV_Counts'][0].path)
        self.assertTrue(chimera_output.exists())
        self.assertGreater(chimera_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
