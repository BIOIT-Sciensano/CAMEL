import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.cutadapt.cutadapt import Cutadapt


class TestCutadapt(CamelTestSuite):
    """
    Tests the cutadapt tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('cutadapt')

    def test_cutadapt_pe(self) -> None:
        """
        Tests basic Cutadapt function in PE mode.
        :return: None
        """
        cutadapt_pe = Cutadapt()
        cutadapt_pe.add_input_files({'FASTQ_PE': [ToolIOFile(TestCutadapt.test_file_dir / 'reads_1.fastq'),
                                                  ToolIOFile(TestCutadapt.test_file_dir / 'reads_2.fastq')]})
        cutadapt_pe.update_parameters(
            minimum_length=20,
            output_basename='output_file',
            quality=10
        )
        cutadapt_pe.run(self.running_dir)
        self.verify_output_files(cutadapt_pe, 'FASTQ_PE', 2)

        self.assertIn('basepair_counts', cutadapt_pe.informs)
        self.assertIn('quality_trimmed', cutadapt_pe.informs['basepair_counts'])
        self.assertIn('quality_trimmed_read1', cutadapt_pe.informs['basepair_counts'])
        self.assertIn('quality_trimmed_read2', cutadapt_pe.informs['basepair_counts'])

    def test_cutadapt_se(self) -> None:
        """
        Tests basic Cutadapt function in SE mode.
        :return: None
        """
        cutadapt_se = Cutadapt()
        cutadapt_se.add_input_files({'FASTQ_SE': [ToolIOFile(TestCutadapt.test_file_dir / 'reads_1.fastq')]})
        cutadapt_se.update_parameters(
            minimum_length=20,
            output_basename='output_file',
            quality=10
        )
        cutadapt_se.run(self.running_dir)
        self.verify_output_files(cutadapt_se, 'FASTQ_SE')

        self.assertIn('basepair_counts', cutadapt_se.informs)
        self.assertIn('quality_trimmed_read1', cutadapt_se.informs['basepair_counts'])

    def test_cutadapt_nextseq(self) -> None:
        """
        Tests Cutadapt NextSeq trimming mode.
        :return: None
        """
        cutadapt_pe = Cutadapt()
        cutadapt_pe.add_input_files({'FASTQ_PE': [
            ToolIOFile(TestCutadapt.test_file_dir / 'reads_1.fastq'),
            ToolIOFile(TestCutadapt.test_file_dir / 'reads_2.fastq')
        ]})

        cutadapt_pe.update_parameters(
            minimum_length=20,
            output_basename='output_file',
            nextseq_trim=10
        )
        cutadapt_pe.run(self.running_dir)
        self.verify_output_files(cutadapt_pe, 'FASTQ_PE', 2)

        self.assertIn('basepair_counts', cutadapt_pe.informs)
        self.assertIn('quality_trimmed', cutadapt_pe.informs['basepair_counts'])
        self.assertIn('quality_trimmed_read1', cutadapt_pe.informs['basepair_counts'])
        self.assertIn('quality_trimmed_read2', cutadapt_pe.informs['basepair_counts'])


if __name__ == '__main__':
    unittest.main()
