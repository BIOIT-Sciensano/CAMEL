import unittest
from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.tools.kma.kmaindex import KMAIndex


class TestKMAIndex(CamelTestSuite):
    """
    Tests the KMA index tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('gene_detection')
    input_fasta = test_file_dir / 'test_input_db.fasta'

    def test_kma_index(self) -> None:
        """
        Tests indexing a FASTA file with KMA.
        :return: None
        """
        kma_index = KMAIndex()
        kma_index.add_input_files({'FASTA': [ToolIOFile(TestKMAIndex.input_fasta)]})
        kma_index.run(self.running_dir)
        self.assertIn('DB', kma_index.tool_outputs, "No DB output generated")

        db_prefix = Path(kma_index.tool_outputs['DB'][0].value)
        index_files = list(db_prefix.parent.glob(f'{db_prefix.name}.*'))
        self.assertGreater(len(index_files), 0, "No KMA index files generated")
        for index_file in index_files:
            self.assertGreater(index_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
