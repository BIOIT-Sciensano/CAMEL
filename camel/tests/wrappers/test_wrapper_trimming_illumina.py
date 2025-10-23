import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.wrappers.trimmingilluminawrapper import TrimmingIlluminaWrapper


class TestWrapperTrimmingIllumina(CamelTestSuite):
    """
    Tests the Illumina trimming wrapper.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_pe = [test_file_dir / 'trimming' / 'reads_1.fastq', test_file_dir / 'trimming' / 'reads_2.fastq']
    fastq_pe_gz = [test_file_dir / 'trimming' / 'reads_1.fastq.gz', test_file_dir / 'trimming' / 'reads_2.fastq.gz']

    def test_trimming_wrapper_illumina_trimmomatic(self) -> None:
        """
        Tests the read trimming wrapper.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingIllumina.fastq_pe)
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)

    def test_trimming_wrapper_illumina_fastp(self) -> None:
        """
        Tests the read trimming wrapper with fastp as trimming method.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingIllumina.fastq_pe, method='fastp')
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)

    def test_trimming_wrapper_illumina_trimmomatic_truseq(self) -> None:
        """
        Tests the read trimming wrapper with TruSeq adapters.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingIllumina.fastq_pe, 'TruSeq3')
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)
        self.assertIn('TruSeq3', wrapper.output.informs_trimming['_command'])

    def test_trimming_wrapper_illumina_trimmomatic_gz(self) -> None:
        """
        Tests the read trimming wrapper with gzipped input.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run(TestWrapperTrimmingIllumina.fastq_pe_gz)
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)


if __name__ == '__main__':
    unittest.main()
