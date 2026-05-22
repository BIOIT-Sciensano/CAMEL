import copy
import tempfile
from pathlib import Path

from Bio.Phylo import NewickIO
from Bio.Phylo.Newick import Tree
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.config import config
from camel.app.loggers import logger
from camel.app.tools.figtree.figtree import FigTree


class NewickUtils:
    """
    This class contains utility functions to work with phylogenetic trees in Newick format.
    """

    @staticmethod
    def calculate_tree_image_height(min_height: int, nb_leaves: int) -> int:
        """
        Calculates the required height for the given tree.
        :param min_height: Minimum height
        :param nb_leaves: Nb. of leaves in the tree
        :return: Image height
        """
        return min_height + 20 * nb_leaves

    @staticmethod
    def create_image_figtree(path_newick: Path, path_config: Path, path_out: Path, width: int, height: int) -> FigTree:
        """
        Creates an image for the given Newick tree using FigTree.
        :param path_newick: Path to Newick file
        :param path_config: Path to FigTree config file
        :param path_out: Output path
        :param width: Image width
        :param height: Image height
        :return: FigTree instance
        """
        logger.info(f"Creating image for tree: {path_newick}")
        with tempfile.TemporaryDirectory(prefix='figtree_', dir=config.dir_temp) as dir_temp:
            figtree = FigTree()
            output_path = Path(dir_temp, 'tree.png')
            figtree.update_parameters(output_path=str(path_out), width=width, height=height)
            path_template = path_config
            figtree.add_input_files({'NWK': [ToolIOFile(path_newick)], 'TXT': [ToolIOFile(path_template)]})
            figtree.run(Path(dir_temp))
        logger.info(f"PNG visualization created: {output_path}")
        return figtree

    @staticmethod
    def count_leaves(newick_in: Path) -> int:
        """
        Counts the number of leaves in a phylogenetic tree.
        :param newick_in: Input Newick file
        :return: Number of leaves
        """
        with open(newick_in) as handle:
            tree = next(NewickIO.parse(handle))
        return tree.count_terminals()

    @staticmethod
    def remove_inner_node_names(tree: Tree) -> Tree:
        """
        Removes inner node names from a tree.
        :param tree: Tree object
        :return: Updated tree
        """
        logger.info("Removing inner node names")
        new_tree = copy.copy(tree)
        for clade in new_tree.find_clades():
            if 'Inner' in clade.name:
                clade.name = ''
        return new_tree

    @staticmethod
    def export_newick_tree(input_tree: Tree, path_out: Path) -> None:
        """
        Saves a tree to a file in Newick format.
        :param input_tree: Input tree
        :param path_out: Output path
        :return: None
        """
        logger.info(f"Saving tree to file: {path_out.name}")
        with path_out.open('w') as handle:
            NewickIO.write([input_tree], handle)
