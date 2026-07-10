import gzip
import re
from pathlib import Path
from typing import Iterable, AbstractSet, Union

from Bio import SeqIO

from camel.app.command.command import Command
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.loggers import logger


class FastqUtils:
    """
    This class contains utility functions to work with FASTQ files.
    """

    @staticmethod
    def count_reads(infile: Path) -> int:
        """
        Count how many reads in a fastq file
        :param infile: file name of the fastq file to count
        :return: number of reads in fastq file
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(infile) else 'cat'
        cmd = f"{cat} {infile} | paste - - - - | wc -l"
        command = Command()
        command.command = cmd
        command.run(infile.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)
        return int(command.stdout.rstrip())

    @staticmethod
    def sort_fastq_by_identifier(infile: Path, outfile: Path, gzip_output: bool = False) -> None:
        """
        Function to sort the reads in a fastq file
        :param infile: file name of the file to sort
        :param outfile: file name of the sorted file
        :param gzip_output: gzip the output file
        :return: None
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(infile) else 'cat'
        gzip_ = '| gzip ' if gzip_output else ''
        cmd = f"{cat} {infile} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" {gzip_}> {outfile}"
        command = Command()
        command.command = cmd
        command.run(infile.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_fastqs_to_interleaved_fastq(
            pe_1_file: Path, pe_2_file: Path, interleaved_file: Path, gzip_output: bool = False) -> None:
        """
        Convert two ordered fastq file into a interleaved fastq file, NO CHECKS are performed. If checks are needed,
        call 'create_paired_end' instead.
        :param pe_1_file: fastq file containing the first group of reads
        :param pe_2_file: fastq file containing the second group of reads
        :param interleaved_file: interleaved fastq file to be generated
        :param gzip_output: gzip the output file
        :return: None
        """
        # (z)cat is used as 'paste - - - - <(zcat filename)' does not work
        cat = 'zcat' if FileSystemHelper.is_gzipped(pe_1_file) else 'cat'
        gzip_ = '| gzip ' if gzip_output else ''
        cmd = f"paste <({cat} {pe_1_file} | paste - - - -) <({cat} {pe_2_file} | paste - - - -) | " \
              f"tr '\\t' '\\n' {gzip_}> {interleaved_file}"
        command = Command()
        command.command = cmd
        command.run(interleaved_file.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    @staticmethod
    def convert_interleaved_fastq_to_individual_fastqs(
            interleaved_file: Path, pe_1_file: Path, pe_2_file: Path, gzip_output: bool = False) -> None:
        """
        Convert an interleaved fastq file into two fastq files each containing one group of reads. Input interleaved
        fastq file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :param gzip_output: gzip the output file
        :return: None
        """
        FastqUtils.split_interleaved_fastq(interleaved_file, pe_1_file, pe_2_file, gzip_output)

    @staticmethod
    def split_interleaved_fastq(
            interleaved_file: Path, pe_1_file: Path, pe_2_file: Path, gzip_output: bool = False, pigz: bool = False) -> None:
        """
        Split an interleaved fastq file into two fastq files each containing one group of reads. Input interleaved fastq
        file must be sorted. No CHECKS are performed.
        :param interleaved_file: interleaved fastq file containing both groups of reads
        :param pe_1_file: fastq file will hold the first group of reads
        :param pe_2_file: fastq file will hold the second group of reads
        :param gzip_output: gzip the output file
        :param pigz: use pigz instead of gzip
        :return: None
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(interleaved_file) else 'cat'
        gzip_ = f"| {'pigz' if pigz else 'gzip'} " if gzip_output else ''
        cmd = f"{cat} {interleaved_file} | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' {gzip_}> {pe_1_file}) |" \
              f" cut -f 5-8 | tr '\t' '\n' {gzip_}> {pe_2_file}"
        command = Command()
        command.command = cmd
        command.run(pe_1_file.resolve().parent)
        if command.stderr != '':
            raise RuntimeError(command.stderr, cmd)

    PATTERN_FQ_PE = r'(.+?)(_S\d+)?(_L\d{3})?[_.]R?1P?(_\d+)?.(fastq|fq)(.gz)?'
    PATTERN_FQ_SE = r'(.+?)(_S\d+)?(_L\d{3})?(_\d+)?.(fastq|fq)(.gz)?'
    PATTERN_FQ_ONT = r'(.+?)(?:_([A-Z0-9-]+)_([A-Z]+v[\d-]+)_(\d+)_(run\d+)_(B\d+))?(?:\.(fastq|fq)(?:\.gz)?)?$'

    @staticmethod
    def get_sample_name(fastq_path: Union[Path, str], pattern: str = PATTERN_FQ_PE) -> str:
        """
        Returns the sample name based on the given reads.
        :param fastq_path: FASTQ path
        :param pattern: Regex to determine the sample name
        :return: Sample name
        """
        basename = FileSystemHelper.make_valid(Path(fastq_path).name)
        m = re.match(pattern, basename, re.IGNORECASE)
        if m:
            return m.group(1)
        raise ValueError(f"Cannot determine sample name from: {basename}")

    @staticmethod
    def get_all_read_names(fastq_path: Path) -> set[str]:
        """
        Retrieves all read names from the given fastq file
        :param fastq_path: Path to the fastq file
        :return: Set with read names
        """
        read_names = set()
        open_fn = gzip.open if FileSystemHelper.is_gzipped(fastq_path) else open
        with open_fn(fastq_path, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                read_names.add(record.id)
        return read_names

    @staticmethod
    @DeprecationWarning
    def _get_read_names(files: Iterable[Path]) -> AbstractSet[str]:
        """
        Returns a set with all the read names in the given files
        :param files: Files for which the read names need to be extracted
        :return: Set of read names
        """
        logger.warning('This function is deprecated')
        read_names = set()
        for file in files:
            read_names = read_names | FastqUtils.get_all_read_names(file)
        return read_names

    @staticmethod
    def count_bases(input_file: Path) -> int:
        """
        Calculates the number of bases in the given input files
        :param input_file: File path
        :return: Number of bases
        """
        cat = 'zcat' if FileSystemHelper.is_gzipped(input_file) else 'cat'
        cmd = f"{cat} {input_file} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c"
        command = Command()
        command.command = cmd
        command.run(Path.cwd())
        return int(command.stdout)

    @staticmethod
    def is_fastq(input_file: Path) -> bool:
        """
        Checks whether the input file is a FASTQ file.
        :param input_file: input file
        :return bool: boolean indicating whether the input file is a FASTQ file or not
        """
        is_gzipped = FileSystemHelper.is_gzipped(input_file)
        open_file = gzip.open(input_file, 'rt') if is_gzipped else open(input_file)
        with open_file as handle:
            fastq = SeqIO.parse(handle, 'fastq')
            try:
                return any(fastq)
            except ValueError:
                return False
