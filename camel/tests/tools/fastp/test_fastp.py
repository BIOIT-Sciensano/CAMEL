import unittest
from pathlib import Path


from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastp.fastp import Fastp
from camel.app.tools.fastp.fastpreporter import FastpReporter


class TestFastp(CamelTestSuite):
    """
    Tests for the fastp tool class.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fastp')
    INPUT_FASTQ_PE = [
        ToolIOFile(test_file_dir / 'reads_illumina_1.fastq.gz'),
        ToolIOFile(test_file_dir / 'reads_illumina_2.fastq.gz')]

    def test_fastp_pe(self) -> None:
        """
        Tests fastp with PE input.
        """
        fastp = Fastp()
        fastp.add_input_files({'FASTQ': TestFastp.INPUT_FASTQ_PE})
        fastp.update_parameters(
            output_name='reads_illumina',
            # Adapter trimming
            detect_adapter_for_pe=True,
            # Leading
            cut_front=True,
            cut_front_window_size=1,
            cut_front_mean_quality=10,
            # Trailing
            cut_tail=True,
            cut_tail_window_size=1,
            cut_tail_mean_quality=10,
            # Sliding window
            cut_right=True,
            cut_right_window_size=4,
            cut_right_mean_quality=20,
            # Minimum length
            length_required=40
        )
        fastp.run(self.running_dir)
        self.verify_output_files(fastp, 'FASTQ_PE', nb_files=2)
        self.verify_output_files(fastp, 'JSON')
        self.verify_output_files(fastp, 'HTML')
        self.assertIn('summary', fastp.informs)

    def test_fastp_se(self) -> None:
        """
        Tests fastp with SE input.
        """
        fastp = Fastp()
        fastp.add_input_files({'FASTQ': [TestFastp.INPUT_FASTQ_PE[0]]})
        fastp.update_parameters(output_name='reads_illumina')
        fastp.run(self.running_dir)
        self.verify_output_files(fastp, 'FASTQ', nb_files=1)
        self.verify_output_files(fastp, 'JSON')
        self.verify_output_files(fastp, 'HTML')
        self.assertIn('summary', fastp.informs)

    def test_fastp_reporter(self) -> None:
        """
        Tests the fastp reporter.
        """
        # Run fastp
        fastp = Fastp()
        fastp.add_input_files({'FASTQ': TestFastp.INPUT_FASTQ_PE})
        fastp.update_parameters(
            output_name='reads_illumina',
            # Adapter trimming
            detect_adapter_for_pe=True,
            # Leading
            cut_front=True,
            cut_front_window_size=1,
            cut_front_mean_quality=10,
            # Trailing
            cut_tail=True,
            cut_tail_window_size=1,
            cut_tail_mean_quality=10,
            # Sliding window
            cut_right=True,
            cut_right_window_size=4,
            cut_right_mean_quality=20,
            # Minimum length
            length_required=40
        )
        fastp.run(self.running_dir)

        # Create dummy FASTQC reports
        path_report_fastqc = Path('fastqc_in.html')
        path_report_fastqc.touch()

        # Create the output report
        reporter = FastpReporter()
        reporter.add_input_files({
            'FASTQ_PE': fastp.tool_outputs['FASTQ_PE'],
            'HTML': fastp.tool_outputs['HTML'],
            'HTML_pre': [ToolIOFile(path_report_fastqc) for _ in range(2)],
            'HTML_post': [ToolIOFile(path_report_fastqc) for _ in range(2)]
        })
        reporter.add_input_informs({'fastp': fastp.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value

        # Save the report in the current directory
        CamelTestSuite.export_report_section(output_section, self.running_dir / 'report')


if __name__ == '__main__':
    unittest.main()
