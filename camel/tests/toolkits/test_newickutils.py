import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.toolkits.phylogeny.newickutils import NewickUtils


class TestNewickUtils(CamelTestSuite):
    """
    Tests the NewickUtils module.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('newick')
    NWK_IN = test_file_dir / 'tree_mst.nwk'
    CONFIG_IN = test_file_dir / 'template_cgmlst_tree.txt'

    def test_create_image_figtree(self) -> None:
        """
        Tests the create_image_figtree function.
        :return: None
        """
        output_png = self.running_dir / 'image.png'
        NewickUtils.create_image_figtree(TestNewickUtils.NWK_IN, TestNewickUtils.CONFIG_IN, output_png, 480, 480)
        self.assertTrue(output_png.exists())
        self.assertGreater(output_png.stat().st_size, 0)

    def test_count_leaves(self) -> None:
        """
        Tests the count_leaves function.
        :return: None
        """
        self.assertEqual(NewickUtils.count_leaves(TestNewickUtils.NWK_IN), 26)


if __name__ == '__main__':
    unittest.main()
