import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.pipelines.klebsiella.kleborate import Kleborate
from camel.app.tools.pipelines.klebsiella.kleboratereporter import KleborateReporter


class TestKleborate(CamelTestSuite):
    """
    Tests the Kleborate tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('kleborate')
    input_fasta = test_file_dir / 'GCF_002248955.1.fna.gz'

    def test_kleborate(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate()
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TSV')

    def test_kleborate_different_confidence(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate()
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(min_kaptive_confidence='Low')
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TSV')

    def test_kleborate_resistance(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate()
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(resistance=True)
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TSV')

    def test_kleborate_all(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate()
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(all=True)
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TSV')

    def test_kleborate_reporter(self) -> None:
        """
        Tests the reporter class for Kleborate.
        :return: None
        """
        # Run Kleborate
        kleborate = Kleborate()
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(all=True)
        kleborate.run(self.running_dir)

        # Create output report
        reporter = KleborateReporter()
        reporter.add_input_files({'TSV': kleborate.tool_outputs['TSV']})
        reporter.add_input_informs({'kleborate': kleborate.informs})
        reporter.run(self.running_dir)

        # Check the output
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
