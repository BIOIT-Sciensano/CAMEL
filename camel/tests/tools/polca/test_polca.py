import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pypolca.pypolca import Pypolca


class TestPypolca(CamelTestSuite):
    """
    Class to test the pypolca tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('pypolca')

    # Create ToolIOFile input files
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_trimmed_1.fastq.gz')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_trimmed_2.fastq.gz')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'consensus.fasta')
    FILE_FASTA_REF_NO_VARIANTS = ToolIOFile(test_file_dir / 'pypolca_corrected.fasta')

    def test_pypolca(self) -> None:
        """
        Tests the Pypolca tool.
        :return: None
        """
        pypolca = Pypolca()
        pypolca.add_input_files({
            'FASTQ_PE': [TestPypolca.FILE_FASTQ_1, TestPypolca.FILE_FASTQ_2],
            'FASTA': [TestPypolca.FILE_FASTA_REF]})
        pypolca.run(self.running_dir)
        self.verify_output_files(pypolca, 'FASTA')

    def test_pypolca_no_variants(self) -> None:
        """
        Tests the Pypolca tool with 0 variants.
        :return: None
        """
        pypolca = Pypolca()
        pypolca.add_input_files({
            'FASTQ_PE': [TestPypolca.FILE_FASTQ_1, TestPypolca.FILE_FASTQ_2],
            'FASTA': [TestPypolca.FILE_FASTA_REF_NO_VARIANTS]})
        pypolca.run(self.running_dir)
        self.verify_output_files(pypolca, 'FASTA')


if __name__ == '__main__':
    unittest.main()
