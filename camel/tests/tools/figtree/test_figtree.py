import unittest
from importlib.resources import files
from pathlib import Path


from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.figtree.figtree import FigTree


class TestFigTree(CamelTestSuite):
    """
    Tests the FigTree tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('newick')
    input_nwk = test_file_dir / 'tree_mst.nwk'

    def test_figtree(self) -> None:
        """
        Tests the FigTree tool.
        :return: None
        """
        path_png_out = self.running_dir / 'tree.png'
        figtree = FigTree()
        figtree.add_input_files({'NWK': [ToolIOFile(TestFigTree.input_nwk)]})
        figtree.update_parameters(output_path=str(path_png_out))
        figtree.run(self.running_dir)
        self.verify_output_files(figtree, 'PNG')

    def test_figtree_with_template(self) -> None:
        """
        Tests the FigTree tool with a template.
        :return: None
        """
        path_png_out = self.running_dir / 'tree.png'
        path_template = Path(str(files('camel').joinpath('resources/tools/figtree/template_cgmlst_tree.txt')))
        figtree = FigTree()
        figtree.add_input_files({'NWK': [ToolIOFile(TestFigTree.input_nwk)], 'TXT': [ToolIOFile(path_template)]})
        figtree.update_parameters(output_path=str(path_png_out))
        figtree.run(self.running_dir)
        self.verify_output_files(figtree, 'PNG')


if __name__ == '__main__':
    unittest.main()
