import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.canu.canu import Canu


class TestCanu(CamelTestSuite):
    """
    Contains tests for the CANU tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('assembly')
    FILE_FQ = ToolIOFile(test_file_dir / 'ont_bsubtilis_small_region.fastq.gz')

    def test_canu(self) -> None:
        """
        Tests Canu 2.2 with ONT data as input.
        """
        canu = Canu()
        canu.add_input_files({'FASTQ': [TestCanu.FILE_FQ]})
        canu.update_parameters(genome_size='10k', minimum_input_coverage=1)
        canu.run(self.running_dir)
        self.verify_output_files(canu, 'FASTA')


if __name__ == '__main__':
    unittest.main()
