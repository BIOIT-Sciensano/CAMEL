import gzip
import re
from pathlib import Path
from typing import Union

from Bio import SeqIO

from camel.app.core.command import Command
from camel.app.core.utils import fileutils


def count_reads(infile: Path) -> int:
    """
    Count how many reads in a fastq file
    :param infile: file name of the fastq file to count
    :return: number of reads in fastq file
    """
    cat = 'zcat' if fileutils.is_gzipped(infile) else 'cat'
    cmd = f"{cat} {infile.name} | paste - - - - | wc -l"
    command = Command(cmd)
    command.run(infile.resolve().parent)
    if command.stderr != '' or command.exit_code != 0:
        raise RuntimeError(command.stderr, cmd)
    return int(command.stdout.rstrip())


PATTERN_FQ_PE = r'(.+?)(_S\d+)?(_L\d{3})?[_.]R?1P?(_\d+)?.(fastq|fq)(.gz)?'
PATTERN_FQ_SE = r'(.+?)(_S\d+)?(_L\d{3})?(_\d+)?.(fastq|fq)(.gz)?'
PATTERN_FQ_ONT = r'(.+?)(?:_([A-Z0-9-]+)_([A-Z]+v[\d-]+)_(\d+)_(run\d+)_(B\d+))?(?:\.(fastq|fq)(?:\.gz)?)?$'

def get_sample_name(fastq_path: Union[Path, str], pattern: str = PATTERN_FQ_PE) -> str:
    """
    Returns the sample name based on the given reads.
    :param fastq_path: FASTQ path
    :param pattern: Regex to determine the sample name
    :return: Sample name
    """
    basename = fileutils.make_valid(Path(fastq_path).name)
    m = re.match(pattern, basename, re.IGNORECASE)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot determine sample name from: {basename}")

def get_all_read_names(fastq_path: Path) -> set[str]:
    """
    Retrieves all read names from the given fastq file
    :param fastq_path: Path to the fastq file
    :return: Set with read names
    """
    read_names = set()
    open_fn = gzip.open if fileutils.is_gzipped(fastq_path) else open
    with open_fn(fastq_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            read_names.add(record.id)
    return read_names

def count_bases(input_file: Path) -> int:
    """
    Calculates the number of bases in the given input files
    :param input_file: File path
    :return: Number of bases
    """
    cat = 'zcat' if fileutils.is_gzipped(input_file) else 'cat'
    cmd = f"{cat} {input_file} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c"
    command = Command(cmd)
    command.run(Path.cwd())
    return int(command.stdout)

def is_fastq(input_file: Path) -> bool:
    """
    Checks whether the input file is a FASTQ file.
    :param input_file: input file
    :return bool: boolean indicating whether the input file is a FASTQ file or not
    """
    is_gzipped = fileutils.is_gzipped(input_file)
    open_file = gzip.open(input_file, 'rt') if is_gzipped else open(input_file)
    with open_file as handle:
        fastq = SeqIO.parse(handle, 'fastq')
        try:
            return any(fastq)
        except ValueError:
            return False
