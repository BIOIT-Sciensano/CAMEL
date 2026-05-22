import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidToolInputError
from camel.app.tools.mothur.mothurmakecontigs import MothurMakeContigs


class TestMakeContigs(CamelTestSuite):
    """
    Tests Mothur make.contigs.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('mothur')
    FILE_FASTQ_FWD = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R1_001.fastq')
    FILE_FASTQ_REV = ToolIOFile(test_file_dir / 'F3D143_S209_L001_R2_001.fastq')

    def test_makecontigs(self) -> None:
        """
        Tests Mothur make.contigs.
        :return: None
        """
        makecontigs = MothurMakeContigs()
        makecontigs.add_input_files({
            'FASTQ_PE': [TestMakeContigs.FILE_FASTQ_FWD, TestMakeContigs.FILE_FASTQ_REV]
        })
        makecontigs.run(self.running_dir)
        self.assertTrue('FASTA_Contig' in makecontigs.tool_outputs, "No contig output generated")
        self.assertTrue('FASTA_Scrap' in makecontigs.tool_outputs, "No scrap output generated")
        self.assertTrue('TSV_Report' in makecontigs.tool_outputs, "No report output generated")
        contig_output = Path(makecontigs.tool_outputs['FASTA_Contig'][0].path)
        self.assertTrue(contig_output.exists())
        self.assertGreater(contig_output.stat().st_size, 0)
        scrap_output = Path(makecontigs.tool_outputs['FASTA_Scrap'][0].path)
        self.assertTrue(scrap_output.exists())
        report_output = Path(makecontigs.tool_outputs['TSV_Report'][0].path)
        self.assertTrue(report_output.exists())
        self.assertGreater(report_output.stat().st_size, 0)

    def test_makecontigs_incorrect_number_of_files(self) -> None:
        """
        Tests whether Mothur make.contigs fails with an incorrect number of input files for FASTQ_PE.
        :return: None
        """
        makecontigs = MothurMakeContigs()
        makecontigs.add_input_files({
            'FASTQ_PE': [TestMakeContigs.FILE_FASTQ_FWD]
        })
        with self.assertRaises(InvalidToolInputError):
            makecontigs.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
