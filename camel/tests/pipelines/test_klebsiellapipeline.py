import unittest

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core import cameltesthelper
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.dbs.dbutils import DBEntry
from camel.app.scriptutils.basescript import basescriptutils
from camel.scripts.klebsiellapipeline import CONFIG_DATA
from camel.scripts.klebsiellapipeline.mainklebsiellapipeline import (
    CUSTOM_ANALYSES, main,
)
from camel.tests import longRunningTest


class TestKlebsiellaPipeline(CamelTestSuite):
    """
    Tests for the Klebsiella pipeline.
    """
    running_dir = None

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe = [
        test_file_dir / 'Kpneumoniae-SRR4046826-ds_1.fastq.gz',
        test_file_dir / 'Kpneumoniae-SRR4046826-ds_2.fastq.gz'
    ]
    input_fastq_se = test_file_dir / 'Kpneumoniae-FAZ63816_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'Kpneumoniae-SRR4046826-ds.fasta'

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
    def test_klebsiella_pipeline_blast(self) -> None:
        """
        Tests the Klebsiella pipeline with blast based detection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestKlebsiellaPipeline.input_fastq_pe[0]),
            str(TestKlebsiellaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'blast',
            '--gene-detection-method', 'blast',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '8'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_kma(self) -> None:
        """
        Tests the Klebsiella pipeline with KMA-based detection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestKlebsiellaPipeline.input_fastq_pe[0]),
            str(TestKlebsiellaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
            '--threads', '8'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_fasta(self) -> None:
        """
        Tests the Klebsiella pipeline using FASTA as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestKlebsiellaPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'blast',
            '--gene-detection-method', 'blast',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a not in ('cgmlst', 'scgmlst')),
            '--threads', '8'
       ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_ont(self) -> None:
        """
        Tests the Klebsiella pipeline using ONT as input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestKlebsiellaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'blast',
            '--gene-detection-method', 'blast',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a not in ('cgmlst', 'scgmlst')),
            '--threads', '8'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_klebsiella_pipeline_ont_kma(self) -> None:
        """
        Tests the Klebsiella pipeline using ONT as input and kma as detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestKlebsiellaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if a not in ('cgmlst', 'scgmlst')),
            '--threads', '8'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
