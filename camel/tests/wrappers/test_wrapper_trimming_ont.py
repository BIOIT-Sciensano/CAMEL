import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.wrappers.trimmingontwrapper import TrimmingONTWrapper


class TestWrapperTrimmingONT(CamelTestSuite):
    """
    Tests the ONT trimming wrapper.
    """
    # Retrieving the test file directory
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_se = test_file_dir / 'trimming' / 'subsampled_ont_13_2500.fastq'
    fastq_se_gz = test_file_dir / 'trimming' / 'subsampled_ont_13_2500.fastq.gz'

    def test_trimming_wrapper_ont(self) -> None:
        """
        Tests the read trimming wrapper.
        :return: None
        """
        wrapper = TrimmingONTWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingONT.fastq_se)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_trimming_wrapper_ont_gz(self) -> None:
        """
        Tests the read trimming wrapper with gzipped input.
        :return: None
        """
        wrapper = TrimmingONTWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingONT.fastq_se_gz)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
