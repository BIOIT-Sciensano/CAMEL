import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.minimap2.minimap2alignment import Minimap2Alignment
from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping


class TestMinimap2(CamelTestSuite):
    """
    Tests the minimap2 tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('minion')
    FILE_REF_GENOME = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'fastq_minion_stec.fastq')
    FILE_FASTA_INFLUENZA = ToolIOFile(
        CamelTestSuite.get_test_file_dir('mash') / 'influenza_a-full_genome.fasta')

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

    def test_minimap2_alignment(self) -> None:
        """
        Tests the minimap2 alignment tool.
        :return: None
        """
        minimap2 = Minimap2Alignment()
        minimap2.add_input_files({
            'FASTA_REF': [TestMinimap2.FILE_FASTA_INFLUENZA],
            'FASTA': [TestMinimap2.FILE_FASTA_INFLUENZA]
        })
        minimap2.run(self.running_dir)
        self.verify_output_files(minimap2, 'PAF')

    def test_minimap2_alignment_cs_tag(self) -> None:
        """
        Tests the minimap2 alignment tool with the cs tag and =/X CIGAR operators enabled.
        :return: None
        """
        minimap2 = Minimap2Alignment()
        minimap2.add_input_files({
            'FASTA_REF': [TestMinimap2.FILE_FASTA_INFLUENZA],
            'FASTA': [TestMinimap2.FILE_FASTA_INFLUENZA]
        })
        minimap2.update_parameters(cs_tag=True, eqx=True, threads=2)
        minimap2.run(self.running_dir)
        self.verify_output_files(minimap2, 'PAF')
        with minimap2.tool_outputs['PAF'][0].path.open() as handle:
            lines = handle.readlines()
        self.assertGreater(len(lines), 0)
        self.assertTrue(all('cs:Z:=' in line for line in lines), 'Missing long cs tag in PAF output')


if __name__ == '__main__':
    unittest.main()
