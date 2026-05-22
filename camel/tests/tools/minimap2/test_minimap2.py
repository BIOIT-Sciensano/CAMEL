import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping


class TestMinimap2(CamelTestSuite):
    """
    Tests the minimap2 tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('minion')
    FILE_REF_GENOME = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'fastq_minion_stec.fastq')

    def test_minimap2(self) -> None:
        """
        Tests the minimap2 tool.
        :return: None
        """
        minimap2 = Minimap2Mapping()
        minimap2.add_input_files({
            'FASTA': [TestMinimap2.FILE_REF_GENOME],
            'FASTQ': [TestMinimap2.FILE_FASTQ]
        })
        minimap2.run(self.running_dir)
        self.assertTrue('SAM' in minimap2.tool_outputs, "No VCF output generated")
        output_file = Path(minimap2.tool_outputs['SAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
