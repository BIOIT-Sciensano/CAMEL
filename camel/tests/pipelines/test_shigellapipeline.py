import unittest

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core import cameltesthelper
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.dbs.dbutils import DBEntry
from camel.app.scriptutils.basescript import basescriptutils
from camel.scripts.shigellapipeline import CONFIG_DATA
from camel.scripts.shigellapipeline.mainshigellapipeline import CUSTOM_ANALYSES, main
from camel.tests import longRunningTest


class TestShigellaPipeline(CamelTestSuite):
    """
    Tests for the Shigella pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe = [
        test_file_dir / 'Shigella-S17BD07654_1.fastq.gz',
        test_file_dir / 'Shigella-S17BD07654_2.fastq.gz'
    ]
    input_fastq_ont = test_file_dir / 'Shigella-SRR29782656_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'Shigella-S17BD07654.fasta'


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
    def test_shigella_pipeline_blast(self) -> None:
        """
        Tests the Shigella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestShigellaPipeline.input_fastq_pe[0]),
            str(TestShigellaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_shigella_pipeline_kma(self) -> None:
        """
        Tests the Shigella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-pe',
            str(TestShigellaPipeline.input_fastq_pe[0]),
            str(TestShigellaPipeline.input_fastq_pe[1]),
            '--input-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_shigella_pipeline_fasta(self) -> None:
        """
        Tests the Shigella pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestShigellaPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_shigella_pipeline_ont_blast(self) -> None:
        """
        Tests the Shigella pipeline with ONT input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestShigellaPipeline.input_fastq_ont),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_shigella_pipeline_ont_kma(self) -> None:
        """
        Tests the Shigella pipeline with ONT input and KMA-based detection with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestShigellaPipeline.input_fastq_ont),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
