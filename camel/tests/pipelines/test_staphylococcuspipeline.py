import unittest

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core import cameltesthelper
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.dbs.dbutils import DBEntry
from camel.app.scriptutils.basescript import basescriptutils
from camel.scripts.staphylococcuspipeline import CONFIG_DATA
from camel.scripts.staphylococcuspipeline.mainstaphylococcuspipeline import (
    CUSTOM_ANALYSES,
    main,
)
from camel.tests import longRunningTest


class TestStaphylococcusPipeline(CamelTestSuite):
    """
    Tests for the Staphylococcus pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Saureus-SRR10393587-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Saureus-SRR10393587-ds_2.fastq.gz'
    ]
    input_fastq_se = test_file_dir / 'pipelines' / 'Saureus-SRR14933399_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'pipelines' / 'Saureus-SRR10393587-ds.fasta'

    def test_dbs(self) -> None:
        """
        Checks if the databases for the pipeline are available.
        :return: None
        """
        data_dbs = cameltesthelper.extract_from_yaml(
            CONFIG_DATA, 'dbs', placeholders={'DB_ROOT': config.dir_db})
        dbs = {key: DBEntry(**data) for key, data in data_dbs.items()}
        self.assertEqual(basescriptutils.check_dbs(dbs), True)

    @longRunningTest()
    def test_illumina_blast(self) -> None:
        """
        Tests the Staphylococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestStaphylococcusPipeline.input_fastq_pe[0]),
            str(TestStaphylococcusPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_kma(self) -> None:
        """
        Tests the Staphylococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestStaphylococcusPipeline.input_fastq_pe[0]),
            str(TestStaphylococcusPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_blast(self) -> None:
        """
        Tests the Staphylococcus pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestStaphylococcusPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_blast(self) -> None:
        """
        Tests the Staphylococcus pipeline with  ONT as input, all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestStaphylococcusPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_kma(self) -> None:
        """
        Tests the Staphylococcus pipeline with ONT input and kma analysis all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestStaphylococcusPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a != 'cgmlst'),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
