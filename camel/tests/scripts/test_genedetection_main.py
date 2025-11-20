import unittest

from camel.app.cli import cliutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.loggers import initialize_logging
from camel.scripts.genedetection.maingenedetection import main


class TestGeneDetection(CamelTestSuite):
    """
    Tests the gene detection workflow.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()

    # FASTA
    input_fasta = test_file_dir / 'workflows' / 'NC_002695.1.fasta'
    input_fasta_galaxy = test_file_dir / 'workflows' / 'dataset_12.dat'

    # FASTQ
    input_fastq_by_key = {
        'illumina': [test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_1.fastq',
                     test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_2.fastq'],
        'ont': [test_file_dir / 'gene_detection' / 'ont' / 'reads_ont.fastq']
    }

    # Special
    input_reads_no_hit = [test_file_dir / 'workflows' / 'ecoli_1.fastq',
                          test_file_dir / 'workflows' / 'ecoli_2.fastq']
    input_reads_raw = [test_file_dir / 'gene_detection' / 'reads-ds_1P.fastq',
                       test_file_dir / 'gene_detection' / 'reads-ds_2P.fastq']
    input_reads_raw_galaxy = [test_file_dir / 'workflows' / 'dataset_fwd_11.dat',
                              test_file_dir / 'workflows' / 'dataset_rev_10.dat']
    input_gene_detection_db = test_file_dir / 'gene_detection' / 'db'

    ###############
    # FASTA input #
    ###############
    def test_gene_detection_blast_fasta(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta),
            '--input-type', 'fasta',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_galaxy(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input in Galaxy format.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--input-type', 'fasta',
            '--fasta-name', TestGeneDetection.input_fasta.name,
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_spaces(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input with spaces in the name.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--fasta-name', '"my reference genome.fasta"',
            '--input-type', 'fasta',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_score(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input with 'score' as blast filtering method.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta',
            str(TestGeneDetection.input_fasta),
            '--input-type', 'fasta',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--blast-min-percent-identity', '95',
            '--blast-min-percent-coverage', '99',
            '--blast-filtering-method', 'score',
            '--blast-score-nb-of-hits', '10',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_overlap(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on FASTA input with 'overlap' as blast filtering method.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta),
            '--input-type', 'fasta',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--blast-min-percent-identity', '95',
            '--blast-min-percent-coverage', '99',
            '--blast-filtering-method', 'overlap',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    ###############
    # FASTQ input #
    ###############
    # BLAST
    def test_gene_detection_blast_illumina(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--input-type', 'illumina',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--assembly-kmers', '33,55',
            '--threads', '4'
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_illumina_fasta_out(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data and saving of the assembly.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        path_fasta_out = self.running_dir / 'report' / 'assembly.fasta'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--input-type', 'illumina',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-fasta', str(path_fasta_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--assembly-kmers', '33,55',
            '--threads', '4'
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)
        self.assertGreater(path_fasta_out.stat().st_size, 0)

    def test_gene_detection_blast_illumina_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--input-type', 'illumina',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads',
            '--assembly-kmers', '33,55',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_ont(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--input-type', 'ont',
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--threads', '4'
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_ont_meta(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on ONT data with meta Flye assembly.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--assembly-flye-meta',
            '--assembly-min-contig-len', '750',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_ont_with_filtering(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on ONT data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--trim-reads',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_ont_reads(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on ONT reads (without assembling step).
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--trim-reads',
            '--blast-reads',
            '--blast-filtering-method', 'score',
            '--blast-score-nb-of-hits', '10',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    # KMA
    def test_gene_detection_kma_illumina(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'illumina',
            '--detection-method', 'kma',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_illumina_trim(self) -> None:
        """
        Tests the gene detection workflow with BLAST detection on Illumina data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'illumina',
            '--detection-method', 'kma',
            '--trim-reads',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_ont(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on ONT data.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--detection-method', 'kma',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_ont_trim(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on ONT data, including trimming.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-se', str(TestGeneDetection.input_fastq_by_key['ont'][0]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'ont',
            '--detection-method', 'kma',
            '--trim-reads',
            '--ont-min-len', '666',
            '--ont-min-qual', '12',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma_apm(self) -> None:
        """
        Tests the gene detection workflow with KMA-based detection on Illumina data, with the apm preset.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe',
            str(TestGeneDetection.input_fastq_by_key['illumina'][0]),
            str(TestGeneDetection.input_fastq_by_key['illumina'][1]),
            '--db', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--input-type', 'illumina',
            '--detection-method', 'kma',
            '--kma-apm', 'p',
        ]
        result = cliutils.invoke(main, args)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
