import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map


class TestBowtie2(CamelTestSuite):
    """
    Tests the Bowtie2 tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('read_mapping')
    path_fasta_in = test_file_dir / 'reference.fasta'
    path_fastq_pe_in = [test_file_dir / 'reads_1.fastq', test_file_dir / 'reads_2.fastq']

    def test_bowtie2_index(self) -> None:
        """
        Tests bowtie2Index.
        :return: None
        """
        bowtie2_index = Bowtie2Index()
        bowtie2_index.add_input_files({'FASTA_REF': [ToolIOFile(TestBowtie2.path_fasta_in)]})
        bowtie2_index.run(self.running_dir)
        self.assertIn('INDEX_GENOME_PREFIX', bowtie2_index.tool_outputs)
        output_file = Path(bowtie2_index.tool_outputs['INDEX_GENOME_PREFIX'][0].value)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_bowtie2_map_pe(self) -> None:
        """
        Tests bowtie2Map in PE mode.
        :return: None
        """
        # Create index
        bowtie2_index = Bowtie2Index()
        bowtie2_index.add_input_files({'FASTA_REF': [ToolIOFile(TestBowtie2.path_fasta_in)]})
        bowtie2_index.run(self.running_dir)

        # Map reads
        bowtie2_map = Bowtie2Map()
        bowtie2_map.add_input_files({
            'FASTQ_PE': [ToolIOFile(TestBowtie2.path_fastq_pe_in[0]), ToolIOFile(TestBowtie2.path_fastq_pe_in[1])],
            'INDEX_GENOME_PREFIX': bowtie2_index.tool_outputs['INDEX_GENOME_PREFIX']
        })
        bowtie2_map.run(self.running_dir)

        # Verify output files
        self.verify_output_files(bowtie2_map, 'SAM')


if __name__ == '__main__':
    unittest.main()
