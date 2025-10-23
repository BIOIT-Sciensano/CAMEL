import logging
import re
import unittest

from Bio import SeqIO
from Bio.Seq import Seq

from camel.app.core.utils import fastautils
from camel.app.core.cameltestsuite import CamelTestSuite


class TestFastautils(CamelTestSuite):
    """
    Tests the FASTA utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')

    def test_rename_sequences_regex(self) -> None:
        """
        Tests the rename_sequences_regex function.
        """
        # Create input FASTA file
        seqs = [SeqIO.SeqRecord(Seq('actg'), id=f'seq_{i}') for i in range(10)]
        path_fasta_in = self.running_dir / 'input.fasta'
        with path_fasta_in.open('w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        logging.debug(f'Dummy FASTA file created: {path_fasta_in}')

        # Replace by regex
        path_out = self.running_dir / 'output.fasta'
        fastautils.rename_sequences_regex(path_fasta_in, path_out, r'seq_(\d+)', r'abc_\g<1>')

        # Check if sequence IDs were replaced
        with path_out.open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                self.assertIsNotNone(re.match(r'abc_\d+', seq.id))

    def test_rename_sequences_with_fasta_file(self) -> None:
        """
        Tests the rename_sequences_with_fasta_file function.
        """
        # Create input FASTA file
        seqs = [SeqIO.SeqRecord(Seq('actg'), id=str(i)) for i in range(10)]
        path_fasta_in = self.running_dir / 'input.fasta'
        with path_fasta_in.open('w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        logging.debug(f'Dummy FASTA file created: {path_fasta_in}')

        # Create reference FASTA file
        seqs = [SeqIO.SeqRecord(Seq('actg'), id=f'seq_{i}') for i in range(10)]
        path_fasta_ref = self.running_dir / 'input.fasta'
        with path_fasta_ref.open('w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        logging.debug(f'Reference FASTA file created: {path_fasta_ref}')

        # Replace sequence ids
        path_out = self.running_dir / 'output.fasta'
        fastautils.rename_sequences_with_fasta_file(path_fasta_in, path_fasta_ref, path_out)

        # Check if sequence IDs were replaced
        with path_out.open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                logging.debug(f'Seq id: {seq.id}')
                self.assertTrue(seq.id.startswith('seq'))

    def test_rename_sequences_with_fasta_file_invalid_input(self) -> None:
        """
        Tests the rename_sequences_with_fasta_file function with invalid input.
        """
        # Create input FASTA file
        seqs_in = [SeqIO.SeqRecord(Seq('actg'), id=str(i)) for i in range(10)]
        path_fasta_in = self.running_dir / 'input.fasta'
        with path_fasta_in.open('w') as handle:
            SeqIO.write(seqs_in, handle, 'fasta')
        logging.debug(f'Dummy FASTA file created: {path_fasta_in}')

        # Create reference FASTA file with a different number of sequences
        seqs_ref = [SeqIO.SeqRecord(Seq('actg'), id=f'seq_{i}') for i in range(4)]
        path_fasta_ref = self.running_dir / 'ref.fasta'
        with path_fasta_ref.open('w') as handle:
            SeqIO.write(seqs_ref, handle, 'fasta')
        logging.debug(f'Reference FASTA file created: {path_fasta_ref}')

        # Replace sequence ids
        path_out = self.running_dir / 'output.fasta'
        with self.assertRaises(ValueError):
            fastautils.rename_sequences_with_fasta_file(path_fasta_in, path_fasta_ref, path_out)

    def test_convert_fasta_to_fastq(self) -> None:
        """
        Test the function to convert a FASTA file to a FASTQ file
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy.fasta'
        output_file = self.running_dir / f"{input_file.stem}.fastq"
        fastautils.convert_fasta_to_fastq(input_file, output_file)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_has_duplicates_false(self) -> None:
        """
        Tests the function that checks whether a FASTA file has duplicate seq IDs.
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy.fasta'
        has_duplicates = fastautils.has_duplicates(input_file)
        self.assertFalse(has_duplicates)

    def test_has_duplicates_true(self) -> None:
        """
        Tests the function that checks whether a FASTA file has duplicate seq IDs.
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy_with_duplicates.fasta'
        has_duplicates = fastautils.has_duplicates(input_file)
        self.assertTrue(has_duplicates)

    def test_count_reads(self) -> None:
        """
        Tests the count_reads function.
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy_with_duplicates.fasta'
        self.assertGreater(fastautils.count_reads(input_file), 0)

    def test_count_bases(self) -> None:
        """
        Tests the count_bases function.
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy_with_duplicates.fasta'
        self.assertGreater(fastautils.count_bases(input_file), 0)

    def test_gc(self) -> None:
        """
        Tests the count_bases function.
        :return: None
        """
        input_file = TestFastautils.test_file_dir / 'toy_with_duplicates.fasta'
        self.assertGreater(fastautils.gc(input_file), 0)


if __name__ == '__main__':
    unittest.main()
