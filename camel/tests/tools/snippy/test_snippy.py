import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.snippy.snippyfasta2bam import SnippyFasta2BAM


class TestSnippy(CamelTestSuite):
    """
    Tests the snippy tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('snippy')
    input_fasta_ref = test_file_dir / 'H37Rv.fasta'
    input_fasta = test_file_dir / 'S32BD03301_contigs_unfilt.fasta'

    def test_snippy_fasta2bam(self) -> None:
        """
        Tests the snippy Fasta2BAM class.
        :return: None
        """
        snippy_fasta2bam = SnippyFasta2BAM()
        snippy_fasta2bam.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSnippy.input_fasta))],
            'FASTA_REF': [ToolIOFile(Path(TestSnippy.input_fasta_ref))]
        })
        snippy_fasta2bam.run(self.running_dir)
        self.verify_output_files(snippy_fasta2bam, 'BAM', 1)


if __name__ == '__main__':
    unittest.main()
