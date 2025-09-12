import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.btyper.btyper import BTyper
from camel.app.tools.btyper.btyperreporter import BTyperReporter


class TestBTyper(CamelTestSuite):
    """
    Tests for BTyper.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('btyper')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bacillus_contigs.fasta')

    def test_btyper(self) -> None:
        """
        Tests BTyper with regular input.
        :return: None
        """
        btyper = BTyper()
        btyper.add_input_files({'FASTA': [TestBTyper.FILE_FASTA]})
        btyper.update_parameters(output_dir=self.running_dir)
        btyper.run(self.running_dir)
        self.verify_output_files(btyper, 'TSV')

    def test_btyper_reporter(self) -> None:
        """
        Tests BTyper followed by the reporter.
        :return: None
        """
        # First run btyper
        btyper = BTyper()
        btyper.add_input_files({'FASTA': [TestBTyper.FILE_FASTA]})
        btyper.update_parameters(output_dir=self.running_dir)
        btyper.run(self.running_dir)
        self.verify_output_files(btyper, 'TSV')

        # Then run the reporter
        reporter = BTyperReporter()
        reporter.add_input_files({'TSV': btyper.tool_outputs['TSV']})
        reporter.add_input_informs({'btyper': btyper.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
