import unittest

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core import cameltesthelper
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.dbs.dbutils import DBEntry
from camel.app.scriptutils.basescript import basescriptutils
from camel.scripts.listeriapipeline import CONFIG_DATA
from camel.scripts.listeriapipeline.mainlisteriapipeline import CUSTOM_ANALYSES, main
from camel.tests import longRunningTest


class TestListeriaPipeline(CamelTestSuite):
    """
    Tests for the Listeria pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Listeria-S16BD02199_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Listeria-S16BD02199_2.fastq.gz'
    ]
    input_fastq_se = test_file_dir / 'pipelines' / 'Listeria-SRR17965220_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'pipelines' / 'Listeria-S16BD02199.fasta'

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
    def test_listeria_pipeline_blast(self) -> None:
        """
        Tests the Listeria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_fasta_out = self.running_dir / 'out' / 'assembly.fasta'
        result = cliutils.invoke(main, [
            '--input-type', 'illumina',
            '--fastq-pe',
            str(TestListeriaPipeline.input_fastq_pe[0]),
            str(TestListeriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-fasta', str(path_fasta_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_kma(self) -> None:
        """
        Tests the Listeria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--input-type', 'illumina',
            '--fastq-pe',
            str(TestListeriaPipeline.input_fastq_pe[0]),
            str(TestListeriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_fasta(self) -> None:
        """
        Tests the Listeria pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fasta', str(TestListeriaPipeline.input_fasta),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst')),
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_ont(self) -> None:
        """
        Tests the Listeria pipeline using FASTA as input with all assays except for cgMLST with ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestListeriaPipeline.input_fastq_se),
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
    def test_listeria_pipeline_kma_ont(self) -> None:
        """
        Tests the Listeria pipeline with all assays except for cgMLST, using KMA with ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--fastq-se', str(TestListeriaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--typing-method', 'kma',
            '--gene-detection-method', 'kma',
            '--analyses', ','.join(a for a in CUSTOM_ANALYSES if not a.startswith('cgmlst'))
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
