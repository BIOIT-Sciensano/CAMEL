import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.polca.polca import Polca


class TestPolca(CamelTestSuite):
    """
    Class to test the polca tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('polca')

    # Create ToolIOFile input files
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_trimmed_1.fastq.gz')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_trimmed_2.fastq.gz')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'consensus.fasta')

    def test_polca(self) -> None:
        """
        Tests the Polca tool.
        :return: None
        """
        polca = Polca()
        polca.add_input_files({
            'FASTQ_PE': [TestPolca.FILE_FASTQ_1, TestPolca.FILE_FASTQ_2],
            'FASTA': [TestPolca.FILE_FASTA_REF]})
        polca.run(self.running_dir)
        self.verify_output_files(polca, 'FASTA')


if __name__ == '__main__':
    unittest.main()
