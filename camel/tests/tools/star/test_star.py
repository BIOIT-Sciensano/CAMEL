import unittest
from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.star.staralign import StarAlign
from camel.app.tools.star.starindex import StarIndex


class TestStar(CamelTestSuite):
    """
    Tests the STAR tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir("star")
    FILE_REF_GENOME = ToolIOFile(test_file_dir / "reference.fasta")
    FILE_FASTQ = ToolIOFile(test_file_dir / "reads.fastq")
    INDEX_DIR = ToolIODirectory(test_file_dir / "GenomeDir")

    def test_star_index(self) -> None:
        """
        Tests StarIndex.
        :return: None
        """
        star_index = StarIndex()
        star_index.add_input_files({
            'FASTA': [TestStar.FILE_REF_GENOME],
        })
        star_index.run(self.running_dir)
        self.assertTrue('INDEX_DIR' in star_index.tool_outputs, "No index generated")
        output_dir = Path(str(star_index.tool_outputs['INDEX_DIR'][0]))
        self.assertTrue(any(output_dir.iterdir()))

    def test_star_align(self) -> None:
        """
        Tests StarAlign (single ends).
        :return: None
        """
        star_align = StarAlign()
        star_align.add_input_files({
            'FASTQ': [TestStar.FILE_FASTQ],
            'INDEX_DIR': [TestStar.INDEX_DIR]
        })
        star_align.run(self.running_dir)
        self.verify_output_files(star_align, 'SAM')


if __name__ == '__main__':
    unittest.main()
