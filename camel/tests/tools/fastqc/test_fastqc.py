import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.fastqc.fastqc import FastQC


class TestFastQC(CamelTestSuite):
    """
    Tests for the FastQC tool class.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fastqc')
    input_fastq_pe_gz = [
        ToolIOFile(test_file_dir / 'reads_illumina_1.fastq.gz'),
        ToolIOFile(test_file_dir / 'reads_illumina_2.fastq.gz')
    ]
    input_fastq_pe = [
        ToolIOFile(test_file_dir / 'reads_illumina_1.fastq'),
        ToolIOFile(test_file_dir / 'reads_illumina_2.fastq'),
    ]

    def test_fastqc(self) -> None:
        """
        Tests FastQC tool with uncompressed input.
        """
        fastqc = FastQC()
        fastqc.add_input_files({'FASTQ': TestFastQC.input_fastq_pe})
        fastqc.run(self.running_dir)
        self.verify_output_files(fastqc, 'HTML', nb_files=2)
        self.verify_output_files(fastqc, 'TXT', nb_files=2)

    def test_fastqc_gz_input(self) -> None:
        """
        Tests FastQC tool with compressed input.
        """
        fastqc = FastQC()
        fastqc.add_input_files({'FASTQ': TestFastQC.input_fastq_pe_gz})
        fastqc.run(self.running_dir)
        self.verify_output_files(fastqc, 'HTML', nb_files=2)
        self.verify_output_files(fastqc, 'TXT', nb_files=2)


if __name__ == '__main__':
    unittest.main()
