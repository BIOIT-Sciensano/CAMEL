import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.integronfinder.integronfinder import IntegronFinder
from camel.app.tools.integronfinder.integronfinderreporter import IntegronFinderReporter


class TestIntegronFinder(CamelTestSuite):
    """
    Tests for the IntegronFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('integron_finder')
    FILE_FASTA = ToolIOFile(test_file_dir / 'GCF_002202175.1.fasta')
    FILE_FASTA_E_COLI = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')

    def test_integron_finder(self) -> None:
        """
        Tests the IntegronFinder tool.
        :return: None
        """
        integron_finder = IntegronFinder()
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

    def test_integron_finder_with_options(self) -> None:
        """
        Tests the IntegronFinder tool.
        :return: None
        """
        integron_finder = IntegronFinder()
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.update_parameters(local_max=True, threads=8)
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

    def test_integron_finder_reporter(self) -> None:
        """
        Tests the IntegronFinder tool with reporter.
        :return: None
        """
        # Run tool
        integron_finder = IntegronFinder()
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

        # Run reporter
        reporter = IntegronFinderReporter()
        reporter.add_input_files({'TSV': integron_finder.tool_outputs['TSV']})
        reporter.add_input_informs({'integron_finder': integron_finder.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
