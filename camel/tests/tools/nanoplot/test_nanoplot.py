import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.nanoplot.nanoplot import NanoPlot


class TestNanoPlot(CamelTestSuite):
    """
    Tests the NanoPlot tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('filtlong')
    input_fastq = test_file_dir / 'reads_nanopore.fastq'
    input_fastq_gz = test_file_dir / 'reads_nanopore.fastq.gz'

    def test_nanoplot(self) -> None:
        """
        Tests the NanoPlot tool.
        :return: None
        """
        nanoplot = NanoPlot()
        nanoplot.add_input_files({'FASTQ': [ToolIOFile(TestNanoPlot.input_fastq)]})
        nanoplot.run(self.running_dir)
        self.verify_output_files(nanoplot, 'TSV')
        self.verify_output_files(nanoplot, 'HTML')
        self.assertIn('number_of_reads', nanoplot.informs)


if __name__ == '__main__':
    unittest.main()
