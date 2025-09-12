from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Seq import Seq

from camel.app.components.blast.blastformat7parser import BlastFormat7Parser

BLASTN_OUTPUT_FORMAT = '"7 pident sseqid sseq slen qseqid qstart qend qseq sstrand"'


@dataclass
class BlastHitStatistics:
    """
    This class contains the statistics for a BLAST hit.
    """
    subject_id: str
    subject_length: int
    subject_sequence: str
    query_id: str
    query_start: int
    query_end: int
    query_sequence: str
    percent_identity: float
    strand: str

    @staticmethod
    def parse_blast_output(output_path: Path) -> list['BlastHitStatistics']:
        """
        Parses a BLAST output file generated with the output format defined above.
        :return: List of parsed BLAST statistics
        """
        return [BlastHitStatistics.create_from_dict(d) for d in BlastFormat7Parser.parse_output_file(output_path)]

    @staticmethod
    def create_from_dict(info: dict[str, Any]) -> 'BlastHitStatistics':
        """
        Creates a BLAST hit statistic from a parsed BLAST output dictionary.
        :param info: Info dictionary
        :return: Blast hit statistics
        """
        try:
            return BlastHitStatistics(
                str(info['sseqid']),
                int(info['slen']),
                info['sseq'],
                str(info['qseqid']),
                int(info['qstart']),
                int(info['qend']),
                str(info['qseq']).replace('-', '') if 'qseq' in info else None,
                float(info['pident']),
                str(info['sstrand']) if 'sstrand' in info else None
            )
        except KeyError as err:
            raise ValueError(f"Key '{err}' missing from blast output")

    @property
    def alignment_length(self) -> int:
        """
        Returns the length of the alignment.
        :return: Length of the alignment
        """
        return len(self.subject_sequence)

    @property
    def subject_coverage(self) -> float:
        """
        Returns the fraction of the subject that is covered by the alignment.
        :return: % subject covered
        """
        return 100.0 * float(self.alignment_length) / self.subject_length

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit
        :return: True if full length hit
        """
        return self.subject_length == self.alignment_length

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect, False otherwise
        """
        return self.is_full_length() and (self.percent_identity == 100.0)

    def is_new_allele(self, min_id: float = 99.0, min_cov: float = 100.0) -> bool:
        """
        Checks if this hit is potentially a novel allele of the locus in the database.
        :param min_id: Min % identity
        :param min_cov: Min % coverage
        :return: True if the allele is new, False otherwise
        """
        if self.is_perfect_hit():
            return False
        if min_cov == 100.0:
            return self.percent_identity >= min_id and self.is_full_length()
        else:
            raise NotImplementedError('Screening for novel alleles with <100% coverage is not implemented')

    def novel_allele_seq(self) -> Seq:
        """
        Returns the sequence of the novel allele.
        :return: Novel allele sequence (if available)
        """
        if not self.is_new_allele():
            raise ValueError('BLAST hit is not a novel allele')
        seq_query = Seq(self.query_sequence)
        return seq_query if self.strand == 'plus' else seq_query.reverse_complement()

    @property
    def gaps(self) -> int:
        """
        Returns the number of gaps in the alignment.
        :return: Number of gaps
        """
        return self.subject_sequence.count('-')

    @property
    def length_statistic(self) -> str:
        """
        Returns the length in the format: {bases_covered}/{subject_length}.
        :return: Length statistic
        """
        return f'{self.alignment_length}/{self.subject_length}'
