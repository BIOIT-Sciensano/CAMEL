import unittest
from pathlib import Path

import yaml


from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
from camel.scripts.stecpipeline import CONFIG_DATA
from camel.scripts.stecpipeline.mainstecpipeline import MainSTECPipeline
from camel.tests import longRunningTest


class TestSTECPipeline(CamelTestSuite):
    """
    Tests for the STEC pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_illumina_pe = [
        test_file_dir / 'STEC-591_S13-ds_1.fastq.gz',
        test_file_dir / 'STEC-591_S13-ds_2.fastq.gz'
    ]
    input_fastq_iontorrent = test_file_dir / 'Ecoli-iontorrent-ERR2019997-ds.fastq.gz'
    input_fastq_ont = test_file_dir / 'STEC-SRR16955601_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'STEC-591_S13-ds.fasta'

    def test_stec_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    def test_stec_pipeline_gene_detection_db(self):
        """
        Checks if the databases for the gene detection are available.
        :return: None
        """
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, db_data in config_data['gene_detection'].items():
            # Check if DB exists
            self.assertGreater(Path(db_data['path']).stat().st_size, 0)

            # Check if metadata and FASTA files can be loaded
            manager = DBManager()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(db_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.tool_outputs), 0)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_stec_pipeline_blast_illumina(self) -> None:
        """
        Tests the STEC pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe',
            str(TestSTECPipeline.input_fastq_illumina_pe[0]),
            str(TestSTECPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_stec_pipeline_kma_illumina(self) -> None:
        """
        Tests the STEC pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe',
            str(TestSTECPipeline.input_fastq_illumina_pe[0]),
            str(TestSTECPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--library', 'TruSeq2'
        ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_stec_pipeline_blast_illumina_with_downsampling(self) -> None:
        """
        Tests the STEC pipeline with Illumina data and downsampling.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe',
            str(TestSTECPipeline.input_fastq_illumina_pe[0]),
            str(TestSTECPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--mlst-warwick',
            '--cov-max', '5.0'
        ]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_stec_pipeline_fasta(self) -> None:
        """
        Tests the STEC pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestSTECPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_stec_pipeline_ont(self) -> None:
        """
        Tests the STEC pipeline using ONT as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestSTECPipeline.input_fastq_ont),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


    @longRunningTest()
    def test_stec_pipeline_kma_ont(self) -> None:
        """
        Tests the STEC pipeline using ONT as input with all assays except for cgMLST and kma as detection method
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
               '--fastq-se', str(TestSTECPipeline.input_fastq_ont),
               '--input-type', 'ont',
               '--output-html', str(path_report_out),
               '--output-dir', str(path_report_out.parent),
               '--output-tsv', str(path_summary_out),
               '--working-dir', str(self.running_dir),
               '--detection-method', 'kma'
           ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
