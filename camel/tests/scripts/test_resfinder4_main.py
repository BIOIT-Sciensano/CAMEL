import unittest
from pathlib import Path

from camel.app.config import config
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.scripts.resfinder.mainresfinder import MainResFinder


class TestResFinder(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('resfinder')
    FILE_FASTA_1 = ToolIOFile(test_file_dir / 'ref_ecoli.fasta')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'salmonella_lt2_ref.fasta')
    FILE_FASTA_3 = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTA_4 = ToolIOFile(test_file_dir / 'enterococcus_unkown_amr_muts.fasta')
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_illumina_2.fastq')
    FILE_FASTA_ENTERO = ToolIOFile(test_file_dir / 'assembly_entero.fasta')
    DB_RESFINDER = Path(config.dir_db, 'resfinder4')

    def test_resfinder_main_fasta(self) -> None:
        """
        Tests the ResFinder main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_3),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--acq-overlap', '40',
            '--point',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Staphylococcus_aureus'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fasta_unknown_phenotypes(self) -> None:
        """
        Tests the ResFinder main script with mutations with unknown phenotypes.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_4),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--min-cov', '60',
            '--threshold', '80',
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fastq(self) -> None:
        """
        Tests the ResFinder main script with FASTQ files.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.FILE_FASTQ_1), str(self.FILE_FASTQ_2),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--point',
            '--acquired',
            '--acq-overlap', '45',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Escherichia_coli'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fastq_kleb(self) -> None:
        """
        Tests the ResFinder main script with FASTQ files for species klebsiella.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.FILE_FASTQ_1), str(self.FILE_FASTQ_2),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--point',
            '--acquired',
            '--acq-overlap', '45',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Klebsiella'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fastq_kleb_only_pointfinder(self) -> None:
        """
        Tests the ResFinder main script (point mutations only) with FASTQ files for species klebsiella.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.FILE_FASTQ_1), str(self.FILE_FASTQ_2),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--point',
            '--acq-overlap', '45',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Klebsiella'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fasta_enterococcus_no_species(self) -> None:
        """
        Tests the ResFinder main script with the Enterococcus input data.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_ENTERO),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--acq-overlap', '40',
            '--min-cov', '60',
            '--threshold', '80',
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
