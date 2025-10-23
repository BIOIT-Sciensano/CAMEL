import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.tools.bwa.bwaindex import BWAIndex
from camel.app.tools.bwa.bwamap import BWAMap


class TestBWA(CamelTestSuite):
    """
    Tests the BWA tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('bwa')

    def test_bwa_index(self) -> None:
        """
        Tests BWAIndex.
        :return: None
        """
        bwa_index = BWAIndex()
        bwa_index.add_input_files({
            'FASTA_REF': [ToolIOFile(TestBWA.test_file_dir / 'reference.fasta')],
        })
        bwa_index.run(self.running_dir)
        self.assertIn('INDEX_GENOME_PREFIX', bwa_index.tool_outputs)
        output_file = Path(bwa_index.tool_outputs['INDEX_GENOME_PREFIX'][0].value)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_bwa_map_pe(self) -> None:
        """
        Tests BWAMap in PE mode.
        :return: None
        """
        bwa_map = BWAMap()
        bwa_map.add_input_files({
            'FASTQ_PE': [
                ToolIOFile(TestBWA.test_file_dir / 'reads_1.fastq'),
                ToolIOFile(TestBWA.test_file_dir / 'reads_2.fastq')
            ],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestBWA.test_file_dir / 'reference.fasta')],
        })
        bwa_map.run(self.running_dir)
        self.assertTrue('SAM' in bwa_map.tool_outputs, "No SAM output generated")
        output_file = Path(bwa_map.tool_outputs['SAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_bwa_map_se(self) -> None:
        """
        Tests BWAMap in SE mode.
        :return: None
        """
        bwa_map = BWAMap()
        bwa_map.add_input_files({
            'FASTQ_SE': [
                ToolIOFile(TestBWA.test_file_dir / 'reads_1.fastq')
            ],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestBWA.test_file_dir / 'reference.fasta')],
        })
        bwa_map.run(self.running_dir)
        self.assertTrue('SAM' in bwa_map.tool_outputs, "No SAM output generated")
        output_file = Path(bwa_map.tool_outputs['SAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
