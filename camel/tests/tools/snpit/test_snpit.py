import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.snpit.snpit import Snpit


class TestSnpit(CamelTestSuite):
    """
    Tests the snpit tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('snpit')
    input_vcf = test_file_dir / 'unfiltered_variants-myco.vcf'

    def test_snpit(self) -> None:
        """
        Tests the snpit tool.
        :return: None
        """
        snpit = Snpit()
        snpit.add_input_files({'VCF': [ToolIOFile(TestSnpit.input_vcf)]})
        snpit.run()
        self.assertIn('species', snpit.informs)


if __name__ == '__main__':
    unittest.main()
