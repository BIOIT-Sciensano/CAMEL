import unittest

from camel.app.core.utils import fastqutils
from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.seqtk.seqtkmergepe import SeqtkMergePE
from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
from camel.app.tools.seqtk.seqtksize import SeqtkSize
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert


class TestSeqtk(CamelTestSuite):
    """
    Initializes the seqtk tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('seqtk')
    FILE_FASTQ_A = ToolIOFile(test_file_dir / 'lm1_1.fastq')
    FILE_FASTQ_B = ToolIOFile(test_file_dir / 'lm1_2.fastq')
    FILE_FASTA = ToolIOFile(test_file_dir / 'lm1_1.fasta')

    def test_seqtk_size(self) -> None:
        """
        Testing SeqKit seq with contigs file.
        :return: None
        """
        seqtk_size = SeqtkSize()
        seqtk_size.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ_A, TestSeqtk.FILE_FASTQ_B]})
        seqtk_size.run(self.running_dir)
        self.assertIn('stats', seqtk_size.informs)
        self.assertEqual(len(seqtk_size.informs['stats']), 2)

    def test_seqtk_convert(self) -> None:
        """
        Testing Seqtk seq -a with fastq file.
        :return: None
        """
        seqtk_convert = SeqtkConvert()
        seqtk_convert.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ_A]})
        seqtk_convert.update_parameters(output_file='sequences.fasta')
        seqtk_convert.run(self.running_dir)
        self.verify_output_files(seqtk_convert, 'FASTA')

    def test_seqtk_seq_fastq(self) -> None:
        """
        Tests seqtk seq with FASTQ input.
        :return: None
        """
        seqtk_seq = SeqtkSeq()
        seqtk_seq.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ_A]})
        seqtk_seq.update_parameters(output_filename='sequences.fastq', sample_fraction='0.5')
        seqtk_seq.run(self.running_dir)
        logger.info(seqtk_seq.tool_outputs)
        self.verify_output_files(seqtk_seq, 'FASTQ')
        self.assertIn('nb_seqs_in', seqtk_seq.informs)
        self.assertIn('nb_seqs_out', seqtk_seq.informs)
        self.assertGreater(seqtk_seq.informs['nb_seqs_in'], seqtk_seq.informs['nb_seqs_out'])

    def test_seqtk_seq_fasta(self) -> None:
        """
        Tests seqtk seq with FASTA input.
        :return: None
        """
        seqtk_seq = SeqtkSeq()
        seqtk_seq.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTA]})
        seqtk_seq.update_parameters(output_filename='sequences.fasta', sample_fraction='0.5')
        seqtk_seq.run(self.running_dir)
        logger.info(seqtk_seq.tool_outputs)
        self.verify_output_files(seqtk_seq, 'FASTA')
        self.assertIn('nb_seqs_in', seqtk_seq.informs)
        self.assertIn('nb_seqs_out', seqtk_seq.informs)
        self.assertGreater(seqtk_seq.informs['nb_seqs_in'], seqtk_seq.informs['nb_seqs_out'])

    def test_seqtk_mergepe(self) -> None:
        """
        Tests seqtk mergepe.
        """
        seqtk_merge_pe = SeqtkMergePE()
        seqtk_merge_pe.add_input_files({'FASTQ_PE': [TestSeqtk.FILE_FASTQ_A, TestSeqtk.FILE_FASTQ_B]})
        seqtk_merge_pe.run(self.running_dir)
        self.verify_output_files(seqtk_merge_pe, 'FASTQ')
        self.assertEqual(
            fastqutils.count_reads(TestSeqtk.FILE_FASTQ_A.path) + fastqutils.count_reads(TestSeqtk.FILE_FASTQ_B.path),
            fastqutils.count_reads(seqtk_merge_pe.tool_outputs['FASTQ'][0].path)
        )


if __name__ == '__main__':
    unittest.main()
