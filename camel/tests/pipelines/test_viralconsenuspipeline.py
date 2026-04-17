import unittest
from pathlib import Path


from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.config import config
from camel.app.core.utils import fastautils
from camel.scripts.viralconsensuspipeline.mainviralconsensuspipeline import main
from camel.tests import longRunningTest


# noinspection DuplicatedCode
class TestViralConsensusPipeline(CamelTestSuite):
    """
    Tests for the viral consensus pipeline.
    """

    dir_testdata = CamelTestSuite.get_test_file_dir('pipelines', 'viral_consensus')
    dir_db = Path(config.dir_db, 'pipelines', 'viral_consensus', 'version_1.1')

    ###################################
    # Illumina - FASTA reference file #
    ###################################
    @longRunningTest()
    def test_illumina_fasta_ref_influenza_a_h1n1(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_a_h1n1_with_json(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_json_out = self.running_dir / 'out' / 'out.json'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--output-json', str(path_json_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_json_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_a_h1n1_with_fasta_export(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        path_fasta_out = self.running_dir / 'out' / 'consensus.fasta'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--output-fasta', str(path_fasta_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(fastautils.count_reads(path_fasta_out), 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_a_h1n1_with_scrubbing(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', 'human_read_scrubbing',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_a_h3n2(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H3N2.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.05-H3N2_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.05-H3N2_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_b_yam(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_b',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_b-YAM.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.08-B_YAM_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.08-B_YAM_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_influenza_b_vic(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_b',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_b-VIC.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.29-B_VIC_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.29-B_VIC_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_sars_cov_2(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'sars_cov_2',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'sars_cov_2-Wuhan-Hu-1.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS2.01_1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS2.01_2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_illumina_fasta_ref_other(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'other',
            '--species-name', 'Respiratory syncytial virus',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'hRSV_A.fasta'),
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'RSV-ERR331022_1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'RSV-ERR331022_2.fastq.gz'),
            '--input-type', 'illumina',
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_fasta_ref_fasta_input(self) -> None:
        """
        Tests the viral consensus pipeline with FASTA input.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fasta', str(TestViralConsensusPipeline.dir_testdata / 'influenza_a-full_genome.fasta'),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_fasta_ref_fasta_input_with_antivirals(self) -> None:
        """
        Tests the viral consensus pipeline with FASTA input and antiviral mutation detection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fasta', str(TestViralConsensusPipeline.dir_testdata / 'influenza_a-h3n2-with_antivirals.fasta'),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    ############################################
    # Illumina - Automatic reference selection #
    ############################################
    @longRunningTest()
    def test_illumina_ref_selection_influenza_a_h1n1(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fastq-pe',
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R1.fastq.gz'),
            str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL2.01-H1N1_R2.fastq.gz'),
            '--input-type', 'illumina',
            '--ref-genome-db', str(TestViralConsensusPipeline.dir_db / 'ref_mash_dbs' / 'influenza_a-gisaid'),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    ##############################
    # ONT - FASTA reference file #
    ##############################
    @longRunningTest()
    def test_ont_fasta_ref_influenza_a_h1n1(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL1.01-H1N1.fastq.gz'),
            '--input-type', 'ont',
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_fasta_ref_influenza_a_h1n1_with_scrubbing(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H1N1.fasta'),
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL1.01-H1N1.fastq.gz'),
            '--input-type', 'ont',
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--analyses', 'human_read_scrubbing',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_fasta_ref_influenza_a_h3n2(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_a',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_a-H3N2.fasta'),
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.INFL1.06_H3N2.fastq.gz'),
            '--input-type', 'ont',
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_ref_selection_influenza_b_missing_segments_fastaref(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data with missing segments.
        Previously the pipeline used to fail on Nextclade when segments were missing but after adding the rule
        'iterative_mapping_add_empty_unselected_segments' in iterative_mapping.smk, the pipeline succeeds.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_b',
            '--fastq-se',
            str(TestViralConsensusPipeline.dir_testdata / 'influenza_b_missing_segments.fastq.gz'),
            '--input-type', 'ont',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'influenza_b-YAM.fasta'),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_fasta_ref_sars_cov_2(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'sars_cov_2',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'sars_cov_2-Wuhan-Hu-1.fasta'),
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS1.01.fastq.gz'),
            '--input-type', 'ont',
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_fasta_ref_sars_cov_2_primer_removal(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'sars_cov_2',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'sars_cov_2-Wuhan-Hu-1.fasta'),
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS1.01.fastq.gz'),
            '--input-type', 'ont',
            '--fasta-primers', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS1.primers.fasta'),
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_fasta_ref_fasta_input_sars_cov_2(self) -> None:
        """
        Tests the viral consensus pipeline with Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'sars_cov_2',
            '--fasta-ref', str(TestViralConsensusPipeline.dir_db / 'ref_genomes' / 'sars_cov_2-Wuhan-Hu-1.fasta'),
            '--fasta', str(TestViralConsensusPipeline.dir_testdata / 'sars_cov_2-BS004897.1.fasta'),
            '--input-type', 'fasta',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    #######################################
    # ONT + automatic reference selection #
    #######################################
    @longRunningTest()
    def test_ont_ref_selection_influenza_b_missing_segments_refsel(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data with missing segments.
        Previously the pipeline used to fail on Nextclade when segments were missing but after adding the rule
        'iterative_mapping_add_empty_unselected_segments' in iterative_mapping.smk, the pipeline succeeds.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'influenza_b',
            '--fastq-se',
            str(TestViralConsensusPipeline.dir_testdata / 'influenza_b_missing_segments.fastq.gz'),
            '--input-type', 'ont',
            '--ref-genome-db', str(TestViralConsensusPipeline.dir_db / 'ref_mash_dbs' / 'influenza_b-gisaid'),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_ont_ref_selection_sars_cov_2(self) -> None:
        """
        Tests the viral consensus pipeline with ONT data & automatic reference selection.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        result = cliutils.invoke(main, [
            '--species', 'sars_cov_2',
            '--fastq-se', str(TestViralConsensusPipeline.dir_testdata / 'ESIB_EQA_2023.SARS1.01.fastq.gz'),
            '--ref-genome-db', str(TestViralConsensusPipeline.dir_db / 'ref_mash_dbs' / 'sars_cov_2-ncbi'),
            '--input-type', 'ont',
            '--clair3-model', str(Path(config.dir_db, 'clair3', 'models', 'ont')),
            '--cov-max', '5000',
            '--cov-max-segment', '500',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
