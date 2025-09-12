import json
import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
from camel.app.tools.pipelines.mycobacterium.amr.amrscreen import AMRScreen


class TestAMRScreen(CamelTestSuite):
    """
    Tests the AMR-screen tool of the Mycobacterium pipeline.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines', 'mycobacterium')
    path_db = Path(config.dir_db, 'pipelines', 'mycobacterium', 'amr', 'who_2023.6')

    def test_amrscreen(self) -> None:
        """
        Tests the AMR-screen tool.
        :return: None
        """
        # Extract variants located in AMR regions
        bcf_filter = BcftoolsFilter()
        bcf_filter.add_input_files({
            'VCF_GZ': [ToolIOFile(TestAMRScreen.test_file_dir / 'variants-S15BD02575-all.vcf.gz')],
            'BED_include': [ToolIOFile(TestAMRScreen.path_db / 'amr_regions.bed')]
        })
        bcf_filter.run(self.running_dir)

        # Run the tool
        amr_screen = AMRScreen()
        amr_screen.add_input_files({
            'VCF': bcf_filter.tool_outputs['VCF'],
            'VCF_filt': bcf_filter.tool_outputs['VCF'],
            'DB': [ToolIODirectory(TestAMRScreen.path_db)],
            'BED': [ToolIOFile(TestAMRScreen.path_db / 'amr_regions.bed')]
        })
        amr_screen.run(self.running_dir)

        # Verify output
        self.verify_output_files(amr_screen, 'TSV')
        self.verify_output_files(amr_screen, 'JSON')
        with amr_screen.tool_outputs['JSON'][0].path.open() as handle:
            results = json.load(handle)
            self.assertGreater(len(results), 0)


if __name__ == '__main__':
    unittest.main()
