import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import ToolExecutionError
from camel.app.tools.itsx.itsx import Itsx


class TestItsx(CamelTestSuite):
    """
    Tests ITSx.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('itsx')
    FILE_FASTA = ToolIOFile(test_file_dir / 'itsx_test.fasta')
    FILE_ERROR = ToolIOFile(test_file_dir / 'empty.fasta')

    def test_itsx(self) -> None:
        """
        Tests ITSx.
        :return: None
        """
        itsx = Itsx()
        itsx.add_input_files({
            'FASTA': [TestItsx.FILE_FASTA]
        })
        itsx.run(self.running_dir)

        for output_key in ['TEXT_Summary', 'TEXT_NoDetection', 'TEXT_Problematic', 'TSV_Positions', 'GRAPH', 'FASTA_Full', 'FASTA_ITS1', 'FASTA_ITS2', 'FASTA_NoDetection']:
            self.assertTrue(output_key in itsx.tool_outputs, f"No {output_key} output generated")
            output_file = Path(itsx.tool_outputs[output_key][0].path)
            self.assertTrue(output_file.exists())
            self.assertGreater(output_file.stat().st_size, 0)

    def test_itsx_outErr(self) -> None:
        """
        Tests the tool execution error raising of the ITSx tool.
        :return: None
        """
        itsx = Itsx()
        itsx.add_input_files({
            'FASTA': [TestItsx.FILE_ERROR]
        })
        with self.assertRaises(ToolExecutionError):
            itsx.run(self.running_dir)

if __name__ == '__main__':
    unittest.main()
