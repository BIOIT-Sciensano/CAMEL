import unittest

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.grapetree.grapetree import GrapeTree


class TestGrapeTree(CamelTestSuite):
    """
    Tests the GrapeTree tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('grapetree')
    input_tsv = test_file_dir / 'allele_matrix.tsv'

    def test_grapetree(self) -> None:
        """
        Tests the GrapeTree tool.
        :return: None
        """
        path_nwk_out = self.running_dir / 'tree.nwk'
        grapetree = GrapeTree()
        grapetree.add_input_files({'TSV': [ToolIOFile(TestGrapeTree.input_tsv)]})
        grapetree.update_parameters(output_path=str(path_nwk_out))
        grapetree.run()
        self.assertIn('NWK', grapetree.tool_outputs)
        self.assertGreater(grapetree.tool_outputs['NWK'][0].size, 0)

    def test_grapetree_nj(self) -> None:
        """
        Tests the GrapeTree tool with the NJ algorithm.
        :return: None
        """
        path_nwk_out = self.running_dir / 'tree.nwk'
        grapetree = GrapeTree()
        grapetree.add_input_files({'TSV': [ToolIOFile(TestGrapeTree.input_tsv)]})
        grapetree.update_parameters(output_path=str(path_nwk_out), method='NJ')
        grapetree.run()
        self.assertIn('NWK', grapetree.tool_outputs)
        self.assertGreater(grapetree.tool_outputs['NWK'][0].size, 0)


if __name__ == '__main__':
    unittest.main()
