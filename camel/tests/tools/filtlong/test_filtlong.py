import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.filtlong.filtlong import Filtlong


class TestFiltlong(CamelTestSuite):
    """
    Tests the Filtlong tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('filtlong')
    input_fastq = test_file_dir / 'reads_nanopore.fastq'
    input_fastq_gz = test_file_dir / 'reads_nanopore.fastq.gz'

    def test_filtlong(self) -> None:
        """
        Tests the Filtlong tool.
        :return: None
        """
        filtlong = Filtlong()
        filtlong.add_input_files({'FASTQ': [ToolIOFile(TestFiltlong.input_fastq)]})
        filtlong.run(self.running_dir)
        self.verify_output_files(filtlong, 'FASTQ')
        self.assertIn('nb_reads_in', filtlong.informs)
        self.assertIn('nb_reads_out', filtlong.informs)

    def test_filtlong_gz_input(self) -> None:
        """
        Tests the Filtlong tool.
        :return: None
        """
        filtlong = Filtlong()
        filtlong.add_input_files({'FASTQ': [ToolIOFile(TestFiltlong.input_fastq_gz)]})
        filtlong.run(self.running_dir)
        self.verify_output_files(filtlong, 'FASTQ')
        self.assertIn('nb_reads_in', filtlong.informs)
        self.assertIn('nb_reads_out', filtlong.informs)

    def test_filtlong_window(self) -> None:
        """
        Tests the Filtlong tool with window quality filtering.
        :return: None
        """
        filtlong = Filtlong()
        filtlong.add_input_files({'FASTQ': [ToolIOFile(TestFiltlong.input_fastq)]})
        filtlong.update_parameters(min_window_q=10)
        filtlong.run(self.running_dir)
        self.verify_output_files(filtlong, 'FASTQ')
        self.assertIn('nb_reads_in', filtlong.informs)
        self.assertIn('nb_reads_out', filtlong.informs)


if __name__ == '__main__':
    unittest.main()
