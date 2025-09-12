import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.genedetection.mainmakegenedetectiondb import MainMakeGeneDetectionDB


class TestMakeGeneDetectionDB(CamelTestSuite):
    """
    Tests the gene detection create DB tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('gene_detection')
    input_fasta = test_file_dir / 'test_input_db.fasta'

    def test_gene_detection_create_db(self) -> None:
        """
        Tests the gene detection create db main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--fasta', str(TestMakeGeneDetectionDB.input_fasta),
            '--working-dir', str(self.running_dir)
        ]
        main = MainMakeGeneDetectionDB(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_gene_detection_create_db_spaces_in_name(self) -> None:
        """
        Tests the gene detection create db main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--fasta', str(TestMakeGeneDetectionDB.input_fasta),
            '--fasta-name', '"spaces in name.fasta"',
            '--working-dir', str(self.running_dir),
            '--threads', '4'
        ]
        main = MainMakeGeneDetectionDB(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
