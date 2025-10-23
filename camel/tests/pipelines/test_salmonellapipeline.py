import unittest
from pathlib import Path

import yaml


from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.scripts.salmonellapipeline import CONFIG_DATA
from camel.scripts.salmonellapipeline.mainsalmonellapipeline import MainSalmonellaPipeline
from camel.tests import longRunningTest


class TestSalmonellaPipeline(CamelTestSuite):
    """
    Tests for the Salmonella pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_illumina_pe = [
        test_file_dir / 'pipelines' / "Salmonella-MB6391-ds_1.fastq.gz",
        test_file_dir / 'pipelines' / "Salmonella-MB6391-ds_2.fastq.gz"
    ]
    input_fastq_se = test_file_dir / 'pipelines' / 'Salmonella-S23BD05337-RBK_ont-ds.fastq.gz'
    input_fasta = test_file_dir / 'pipelines' / "Salmonella-MB6391-ds.fasta"

    def test_salmonella_pipeline_typing_db(self) -> None:
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

    def test_salmonella_pipeline_gene_detection_db(self):
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
    def test_salmonella_pipeline_blast_illumina(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_fasta_out = self.running_dir / 'out' / 'assembly.fasta'
        args = [
            '--fastq-pe', str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-fasta', str(path_fasta_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    @longRunningTest()
    def test_salmonella_pipeline_blast_illumina_with_json(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_json_out = self.running_dir / 'out' / 'output.json'
        args = [
            '--fastq-pe', str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-json', str(path_json_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_json_out.stat().st_size, 0)

    @longRunningTest()
    def test_salmonella_pipeline_kma_illumina(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestSalmonellaPipeline.input_fastq_illumina_pe[0]),
            str(TestSalmonellaPipeline.input_fastq_illumina_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
        ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_salmonella_pipeline_fasta(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestSalmonellaPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_salmonella_pipeline_ont(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST , ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestSalmonellaPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_salmonella_pipeline_kma_ont(self) -> None:
        """
        Tests the Salmonella pipeline with all assays except for cgMLST, using KMA and ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-se', str(TestSalmonellaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
        ] + [f"--{a.replace('_', '-')}" for a in MainSalmonellaPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainSalmonellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
