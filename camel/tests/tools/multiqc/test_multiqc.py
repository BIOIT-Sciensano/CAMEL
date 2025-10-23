import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import ToolExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.multiqc.multiqc import MultiQC


class TestMultiQC(CamelTestSuite):
    """
    Tests the bcftools tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('multiqc')

    def test_multiqc(self) -> None:
        """
        Tests basic multiqc function.
        :return: None
        """
        multiqc = MultiQC()
        multiqc.add_input_files({
            'TXT_log': [
                ToolIOFile(TestMultiQC.test_file_dir / "GC088815.agg.alignment_summary_metrics"),
                ToolIOFile(TestMultiQC.test_file_dir / "GC088816.agg.alignment_summary_metrics"),
                ToolIOFile(TestMultiQC.test_file_dir / "GC088817.agg.alignment_summary_metrics")]
        })
        multiqc.update_parameters(output_directory=self.running_dir)
        multiqc.run(self.running_dir)
        self.verify_output_files(multiqc, 'HTML')

        self.assertIn('report_name', multiqc.informs)
        self.assertIn('data_dir', multiqc.informs)

    def test_multiqc_error(self) -> None:
        """
        Tests multiqc with broken custom file.
        :return: None
        """
        multiqc_error = MultiQC()
        multiqc_error.add_input_files({
            'TXT_log': [ToolIOFile(TestMultiQC.test_file_dir / "QC_summary_mqc.yml")]
        })
        multiqc_error.update_parameters(output_directory=self.running_dir)

        with self.assertRaises(ToolExecutionError):
            multiqc_error.run(self.running_dir)

    def test_multiqc_ignore(self) -> None:
        """
        Tests MultiQC with ignore parameter
        :return: None
        """
        multiqc = MultiQC()
        multiqc.add_input_files({
            'TXT_log': [ToolIOFile(TestMultiQC.test_file_dir / "GC088815.agg.alignment_summary_metrics")]
        })
        multiqc.update_parameters(
            output_directory=self.running_dir,
            ignore=f"{str(TestMultiQC.test_file_dir)}/GC088816.agg.alignment_summary_metrics;{str(TestMultiQC.test_file_dir)}/GC088817.agg.alignment_summary_metrics"
        )
        multiqc.run(self.running_dir)
        self.verify_output_files(multiqc, 'HTML')

        self.assertIn('report_name', multiqc.informs)
        self.assertIn('data_dir', multiqc.informs)

if __name__ == '__main__':
    unittest.main()
