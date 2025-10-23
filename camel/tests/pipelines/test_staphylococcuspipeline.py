import unittest
from pathlib import Path

import yaml


from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.scripts.staphylococcuspipeline import CONFIG_DATA
from camel.scripts.staphylococcuspipeline.mainstaphylococcuspipeline import MainStaphylococcusPipeline
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

    def test_staphylococcus_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing']['dbs'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = TypingDBLoader()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    def test_staphylococcus_pipeline_gene_detection_db(self) -> None:
        """
        Checks if the databases for the gene detection are available.
        :return: None
        """
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, db_data in config_data['gene_detection']['dbs'].items():
            # Check if scheme exists
            self.assertGreater(Path(db_data['path']).stat().st_size, 0)

            # Check if metadata and FASTA files can be loaded
            manager = DBManager()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(db_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.tool_outputs), 0)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_staphylococcus_pipeline_blast(self) -> None:
        """
        Tests the Staphylococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestStaphylococcusPipeline.input_fastq_pe[0]),
            str(TestStaphylococcusPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainStaphylococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainStaphylococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_staphylococcus_pipeline_kma(self) -> None:
        """
        Tests the Staphylococcus pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestStaphylococcusPipeline.input_fastq_pe[0]), str(TestStaphylococcusPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainStaphylococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainStaphylococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_staphylococcus_pipeline_fasta(self) -> None:
        """
        Tests the Staphylococcus pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestStaphylococcusPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainStaphylococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainStaphylococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_staphylococcus_pipeline_ont(self) -> None:
        """
        Tests the Staphylococcus pipeline with  ONT as input, all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestStaphylococcusPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainStaphylococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainStaphylococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_staphylococcus_pipeline_kma_ont(self) -> None:
        """
        Tests the Staphylococcus pipeline with ONT input and kma analysis all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-se', str(TestStaphylococcusPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainStaphylococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainStaphylococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
