import unittest
from pathlib import Path

import yaml


from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.scripts.listeriapipeline import CONFIG_DATA
from camel.scripts.listeriapipeline.mainlisteriapipeline import MainListeriaPipeline
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

    def test_listeria_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
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

    def test_listeria_gene_detection_db(self):
        """
        Checks if the databases for the gene detection are available.
        :return: None
        """
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, db_data in config_data['gene_detection'].items():
            # Check if scheme exists
            self.assertGreater(Path(db_data['path']).stat().st_size, 0)

            # Check if metadata and FASTA files can be loaded
            manager = DBManager()
            manager.add_input_files({'DIR': [ToolIODirectory(Path(db_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.tool_outputs), 0)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_listeria_pipeline_blast(self) -> None:
        """
        Tests the Listeria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_fasta_out = self.running_dir / 'out' / 'assembly.fasta'
        args = [
            '--fastq-pe', str(TestListeriaPipeline.input_fastq_pe[0]),
            str(TestListeriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-fasta', str(path_fasta_out),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainListeriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainListeriaPipeline(args)
        main.run()
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
        args = [
            '--fastq-pe', str(TestListeriaPipeline.input_fastq_pe[0]),
            str(TestListeriaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainListeriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainListeriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_fasta(self) -> None:
        """
        Tests the Listeria pipeline using FASTA as input with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fasta', str(TestListeriaPipeline.input_fasta),
                   '--input-type', 'fasta',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainListeriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainListeriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_ont(self) -> None:
        """
        Tests the Listeria pipeline using FASTA as input with all assays except for cgMLST with ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
                   '--fastq-se', str(TestListeriaPipeline.input_fastq_se),
                   '--input-type', 'ont',
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir)
               ] + [f"--{a.replace('_', '-')}" for a in MainListeriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainListeriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_listeria_pipeline_kma_ont(self) -> None:
        """
        Tests the Listeria pipeline with all assays except for cgMLST, using KMA with ONT input
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-se', str(TestListeriaPipeline.input_fastq_se),
            '--input-type', 'ont',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainListeriaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainListeriaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
