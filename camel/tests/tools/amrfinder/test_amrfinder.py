import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.amrfinder.amrfinder import AMRFinder
from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter


class TestAMRFinder(CamelTestSuite):
    """
    Tests the AMRFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('amrfinder')
    FILE_FASTA = ToolIOFile(test_file_dir / 'test_dna.fa')
    DIR_DB = config.dir_db / 'amrfinder' / 'v4' / 'latest'

    def test_amrfinder(self) -> None:
        """
        Tests the AMRFinder tool with default options.
        :return: None
        """
        amrfinder = AMRFinder()
        amrfinder.add_input_files({'FASTA': [TestAMRFinder.FILE_FASTA], 'DIR': [ToolIODirectory(TestAMRFinder.DIR_DB)]})
        amrfinder.run(self.running_dir)
        self.verify_output_files(amrfinder, 'TSV')

    def test_amrfinder_params(self) -> None:
        """
        Tests the AMRFinder tool with parameters.
        :return: None
        """
        amrfinder = AMRFinder()
        amrfinder.add_input_files({'FASTA': [TestAMRFinder.FILE_FASTA], 'DIR': [ToolIODirectory(TestAMRFinder.DIR_DB)]})
        amrfinder.update_parameters(min_cov=0.9, min_ident=0.9)
        amrfinder.run(self.running_dir)
        self.verify_output_files(amrfinder, 'TSV')

    def test_amrfinder_escherichia(self) -> None:
        """
        Tests the AMRFinder tool with the 'organisms' parameter.
        :return: None
        """
        amrfinder = AMRFinder()
        amrfinder.add_input_files({'FASTA': [TestAMRFinder.FILE_FASTA], 'DIR': [ToolIODirectory(TestAMRFinder.DIR_DB)]})
        amrfinder.update_parameters(organism='Escherichia')
        amrfinder.run(self.running_dir)
        self.verify_output_files(amrfinder, 'TSV')

    def test_amrfinder_reporter(self) -> None:
        """
        Tests the AMRFinder reporter class.
        :return: None
        """
        # Run AMRFinder
        amrfinder = AMRFinder()
        amrfinder.add_input_files({'FASTA': [TestAMRFinder.FILE_FASTA], 'DIR': [ToolIODirectory(TestAMRFinder.DIR_DB)]})
        amrfinder.update_parameters(organism='Escherichia')
        amrfinder.run(self.running_dir)
        self.verify_output_files(amrfinder, 'TSV')

        # Run the AMRFinder reporter
        reporter = AMRFinderReporter()
        reporter.add_input_files({'TSV': amrfinder.tool_outputs['TSV']})
        reporter.add_input_informs({'amrfinder': amrfinder.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
