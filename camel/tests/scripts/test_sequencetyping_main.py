import unittest
from pathlib import Path

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.scripts.sequencetyping.mainsequencetyping import MainSequenceTyping


class TestSequenceTyping(CamelTestSuite):
    """
    Tests the sequence typing tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('typing')
    input_db_nucl = test_file_dir / 'scheme_mlst_neisseria'
    input_db_protein = test_file_dir / 'scheme_pora_neisseria'
    input_db_mixed = test_file_dir / 'scheme_fhbp_neisseria'
    input_fasta = test_file_dir / 'neisseria_mc58.fasta'
    input_fasta_new_allele = test_file_dir / 'neisseria_mc58_new_allele.fasta'
    input_fasta_multi_perfect = test_file_dir / 'neisseria_multi_perfect.fasta'
    input_typing_reads = {
        'illumina': [test_file_dir / 'S15BD05018_S58_L001_1.fastq', test_file_dir / 'S15BD05018_S58_L001_2.fastq'],
        'ont': [test_file_dir / 'Z1269_RPB-subset.fastq.gz']
    }

    def test_typing_fasta_blast_nucl(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_fasta_blast_nucl_tsv_out(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions) and the tabular output
        option.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        output_file_tsv = Path(self.running_dir) / 'report' / 'typing_out.tsv'
        args = [
            '--fasta', str(self.input_fasta),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--output-tsv', str(output_file_tsv),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)
        self.assertGreater(output_file_tsv.stat().st_size, 0)

    def test_typing_fasta_blast_nucl_new_allele(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        The input file was modified to include:
        - A novel allele with a single SNP
        - A novel allele that is also a multi-hit
        - A missing locus from the MLST scheme
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta_new_allele),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_fasta_blast_nucl_multi_perfect_hit(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        The input file was modified to include:
        - A novel allele with a single SNP
        - A novel allele that is also a multi-hit
        - A missing locus from the MLST scheme
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta_multi_perfect),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_fasta_blast_peptide(self) -> None:
        """
        Tests sequence typing using BLAST with a peptide scheme.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_protein),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_fasta_blast_mixed(self) -> None:
        """
        Tests sequence typing using BLAST with a mixed scheme (DNA & peptide loci).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_illumina_kma_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
            '--input-type', 'illumina',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_illumina_kma_mixed(self) -> None:
        """
        Tests sequence typing using KMA with a mixed scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
            '--input-type', 'illumina',
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_ont_blast_nucl(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        Includes de novo assembly
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['ont'][0]),
            '--input-type', 'ont',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_ont_kma_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['ont'][0]),
            '--input-type', 'ont',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_ont_kma_nucl_trim(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['ont'][0]),
            '--input-type', 'ont',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8',
            '--trim-reads'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_fasta_mist_nucl(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--input-type', 'fasta',
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'mist',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
