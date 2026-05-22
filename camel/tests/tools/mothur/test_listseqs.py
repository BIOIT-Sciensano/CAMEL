import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.mothur.mothurlistseqs import MothurListSeqs


class TestListSeqs(CamelTestSuite):
    """
    Tests Mothur list.seqs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_TSV_LIST = ToolIOFile(test_file_dir / 'ITS.list')

    def test_listseqs(self) -> None:
        """
        Tests Mothur list.seqs.
        :return: None
        """
        listseqs = MothurListSeqs()
        listseqs.add_input_files({
            'TSV_List': [TestListSeqs.FILE_TSV_LIST]
        })
        listseqs.run(self.running_dir)
        self.assertTrue('TSV_Accnos' in listseqs.tool_outputs, "No accnos file generated")
        accnos_output = Path(listseqs.tool_outputs['TSV_Accnos'][0].path)
        self.assertTrue(accnos_output.exists())
        self.assertGreater(accnos_output.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
