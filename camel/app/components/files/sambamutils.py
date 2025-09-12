from pathlib import Path

import pysam
from Bio import SeqIO

from camel.app.loggers import logger


class SAMBAMutils:
    """
    Helper to perform SAM BAM file related functions
    """

    @staticmethod
    def get_record_count(path_in: Path) -> int:
        """
        Count the total number of records in a SAM/BAM file
        :param path_in: Input SAM/BAM file
        :return: the number of records in SAM/BAM file
        """
        mode = 'rb' if path_in.name.endswith('.bam') else 'r'
        with pysam.AlignmentFile(str(path_in), mode=mode) as bam:
            return bam.count(until_eof=True)

    @staticmethod
    def is_empty(path_in: Path) -> bool:
        """
        Check whether a SAM/BAM file contains only header
        :param path_in: Input SAM/BAM file
        :return: True if infile contains no read records
        """
        return SAMBAMutils.get_record_count(path_in) == 0

    @staticmethod
    def create_empty(path_out: Path, fasta_in: Path, compress=False) -> None:
        """
        Creates an empty BAM/SAM file based on the input reference genome.
        :param path_out: Output path
        :param fasta_in: Reference FASTA file
        :param compress: If true, output is generated in BAM format, otherwise in SAM format
        :return: None
        """
        with fasta_in.open() as handle:
            sq_list = [{'SN': seq.id, 'LN': len(seq)} for seq in SeqIO.parse(handle, 'fasta')]
        mode = 'wh' if not compress else 'wb'
        with pysam.AlignmentFile(str(path_out), mode=mode, header={'SQ': sq_list}) as _:
            pass
        logger.info(f'Empty BAM/SAM file created: {path_out}')
