import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.viral_consensus.callmultiallelicsites import CallMultiAllelicSites


class TestCallMultiAllelicSites(CamelTestSuite):
    """
    Tests for the multi-allelic site calling.
    """

    path_vcf_in = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus') / 'pileup_multi_allelic.bcf'

    def test_call_multi_allelic_sites(self) -> None:
        """
        Tests the call multi-allelic site calling.
        :return: None
        """
        call_multi_allelic_sites = CallMultiAllelicSites()
        call_multi_allelic_sites.add_input_files({'VCF': [ToolIOFile(TestCallMultiAllelicSites.path_vcf_in)]})
        call_multi_allelic_sites.run(self.running_dir)
        self.assertIn('nb_sites', call_multi_allelic_sites.informs)
        self.verify_output_files(call_multi_allelic_sites, 'TSV')


if __name__ == '__main__':
    unittest.main()
