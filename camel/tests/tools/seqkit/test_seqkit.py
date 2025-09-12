import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqkit.seqkitseq import SeqkitSeq
from camel.app.tools.seqkit.seqkitstats import SeqkitStats
from camel.app.tools.seqkit.seqkitsplit2 import SeqkitSplit2


class TestSeqkit(CamelTestSuite):
    """
    Tests for the seqkit module.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('seqkit')
    FILE_FASTA = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')

    def test_seqkit_seq_fasta(self) -> None:
        """
        Testing SeqKit seq with contigs file.
        """
        seqkit_seq = SeqkitSeq()
        seqkit_seq.add_input_files({'FASTA': [TestSeqkit.FILE_FASTA]})
        seqkit_seq.update_parameters(output_filename='sequences.fasta', max_length=100000)
        seqkit_seq.run(self.running_dir)
        self.verify_output_files(seqkit_seq, 'FASTA')
        self.assertIn('nb_seqs_in', seqkit_seq.informs)
        self.assertIn('nb_seqs_out', seqkit_seq.informs)
        self.assertGreater(seqkit_seq.informs['nb_seqs_in'], seqkit_seq.informs['nb_seqs_out'])

    def test_seqkit_seq_fastq(self) -> None:
        """
        Testing SeqKit seq with paired-end fastq reads.
        """
        seqkit_seq = SeqkitSeq()
        seqkit_seq.add_input_files({'FASTQ': [TestSeqkit.FILE_FASTQ]})
        seqkit_seq.update_parameters(max_length=150, min_qual=7)
        seqkit_seq.run(self.running_dir)
        self.verify_output_files(seqkit_seq, 'FASTQ')
        self.assertIn('nb_seqs_in', seqkit_seq.informs)
        self.assertIn('nb_seqs_out', seqkit_seq.informs)
        self.assertGreater(seqkit_seq.informs['nb_seqs_in'], seqkit_seq.informs['nb_seqs_out'])

    def test_seqkit_stats_fastq(self) -> None:
        """
        Testing SeqKit stats with single-end fastq reads.
        """
        seqkit_stats = SeqkitStats()
        seqkit_stats.add_input_files({'FASTQ': [TestSeqkit.FILE_FASTQ]})
        seqkit_stats.run(self.running_dir)
        self.verify_output_files(seqkit_stats, 'TSV')

    def test_seqkit_split2_fasta(self) -> None:
        """
        Testing SeqKit split2 with fasta file.
        """
        seqkit_split = SeqkitSplit2()
        seqkit_split.add_input_files({'FASTA': [TestSeqkit.FILE_FASTA]})
        seqkit_split.update_parameters(by_part=5, output_dir=self.running_dir / 'parts')
        seqkit_split.run(self.running_dir)
        self.verify_output_files(seqkit_split, 'FASTA', nb_files=5)


if __name__ == '__main__':
    unittest.main()
