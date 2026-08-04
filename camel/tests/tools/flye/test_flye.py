import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.flye.flye import Flye


class TestFlye(CamelTestSuite):
    """
    Contains tests for Flye.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('assembly')
    FILE_FQ = ToolIOFile(test_file_dir / 'ont_bsubtilis_small_region.fastq.gz')

    def test_flye(self) -> None:
        """
        Tests Flye 2.9.4 with ONT data as input.
        """
        flye = Flye()
        flye.add_input_files({'FASTQ': [TestFlye.FILE_FQ]})
        flye.update_parameters(nano_corr=False, nano_hq='', genome_size='10k')
        flye.run(self.running_dir)
        self.verify_output_files(flye, 'FASTA')


if __name__ == '__main__':
    unittest.main()
