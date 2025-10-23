import unittest

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.tools.polypolish.polypolish import Polypolish
from camel.app.tools.polypolish.polypolishinsertfilter import PolypolishInsertFilter


class TestPolypolish(CamelTestSuite):
    """
    Initializes the Polypolish testing tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('polypolish')

    # Create ToolIOFile input files
    FILE_SAM = ToolIOFile(test_file_dir / 'S1_ilmn.sam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'assembly.fasta')

    # Create ToolIOFile for multiple SAM inputs
    FILE_SAM_1 = ToolIOFile(test_file_dir / 'alignment_1_filtered.sam')
    FILE_SAM_2 = ToolIOFile(test_file_dir / 'alignment_2_filtered.sam')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'input_genome.fasta')

    # Create ToolIOFile for insertm filter
    UNFILTERED_SAM_1 = ToolIOFile(test_file_dir / 'alignment_1.sam')
    UNFILTERED_SAM_2 = ToolIOFile(test_file_dir / 'alignment_2.sam')

    def test_polypolish(self) -> None:
        """
        Tests Polypolish.
        :return: None
        """
        polypolish = Polypolish()
        polypolish.add_input_files({
            'SAM': [TestPolypolish.FILE_SAM],
            'FASTA': [TestPolypolish.FILE_FASTA_REF]})
        polypolish.run(self.running_dir)
        self.verify_output_files(polypolish, 'FASTA')

    def test_polypolish_multiple_sam(self) -> None:
        """
        Tests Polypolish.
        :return: None
        """
        polypolish = Polypolish()
        polypolish.add_input_files({
            'SAM': [TestPolypolish.FILE_SAM_1, TestPolypolish.FILE_SAM_2],
            'FASTA': [TestPolypolish.FILE_FASTA_2]})
        polypolish.run(self.running_dir)
        self.verify_output_files(polypolish, 'FASTA')

    def test_polypolish_insert_filter(self) -> None:
        """
        Tests PolypolishInsertFilter.
        :return: None
        """
        polypolish_insert_filter = PolypolishInsertFilter()
        polypolish_insert_filter.add_input_files({
            'SAM': [TestPolypolish.UNFILTERED_SAM_1, TestPolypolish.UNFILTERED_SAM_2]})
        polypolish_insert_filter.run(self.running_dir)
        self.verify_output_files(polypolish_insert_filter, 'SAM', nb_files=2)


if __name__ == '__main__':
    unittest.main()
