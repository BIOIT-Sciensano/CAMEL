import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.scripts.confindr.mainconfindr import main


class TestConFindr(CamelTestSuite):
    """
    Tests the ConFindr tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('confindr')
    input_pe_reads = [
        test_file_dir / 'illumina_enterococcus_1.fastq.gz',
        test_file_dir / 'illumina_enterococcus_2.fastq.gz'
    ]
    input_se_reads = test_file_dir / 'minion_reads-ecoli.fastq'
    db = config.dir_db / 'confindr' / '0.8.1'

    def test_confindr_main_script_se(self) -> None:
        """
        Tests the ConFinder main script.
        :return: None
        """
        dir_out = self.running_dir / 'out'
        dir_out.mkdir()
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestConFindr.input_se_reads),
            '--db', str(TestConFindr.db),
            '--working-dir', str(self.running_dir),
            '--output-html', str(dir_out / 'report.html'),
            '--output-dir', str(dir_out),
            '--quality-cutoff', '15',
            '--base-cutoff', '5',
            '--base-percentage-cutoff', '10',
            '--min-matching-hashes', '200',
            '--input-type', 'ont',
            '--rmlst'
        ])
        self.assertEqual(result.exit_code, 0)

    def test_confindr_main_script_pe(self) -> None:
        """
        Tests the ConFinder main script.
        :return: None
        """
        dir_out = self.running_dir / 'out'
        dir_out.mkdir()
        result = cliutils.invoke(main, [
            '--fastq-pe', str(TestConFindr.input_pe_reads[0]), str(TestConFindr.input_pe_reads[1]),
            '--fastq-pe-names', TestConFindr.input_pe_reads[0].name, TestConFindr.input_pe_reads[1].name,
            '--db', str(TestConFindr.db),
            '--working-dir', str(self.running_dir),
            '--output-html', str(dir_out / 'report.html'),
            '--output-dir', str(dir_out),
            '--quality-cutoff', '15',
            '--base-cutoff', '5',
            '--base-percentage-cutoff', '10',
            '--min-matching-hashes', '200',
            '--input-type', 'illumina',
            '--rmlst'
        ])
        self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()
