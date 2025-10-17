import re
from pathlib import Path
from typing import Union, Optional

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from camel.app.core.command import Command


def read_as_dict(fasta: Path) -> dict:
    """
    Read in fasta file as a dictionary keyed by sequence id. More efficient for small files.
    :param fasta: fasta file to read
    :return: a sequence record dictionary
    """
    with fasta.open() as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

def write(sequences: list[Union[SeqIO.SeqRecord, str]], fasta: Path) -> None:
    """
    Write a list of sequence records into fasta file
    :param sequences: list of sequences as SeqRecord object
    :param fasta: fasta file to write into
    :return: None
    """
    with open(fasta, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def count_reads(infile: Path) -> int:
    """
    Counts the number of reads in a FASTA file.
    :param infile: Input FASTA file
    :return: Number of reads
    """
    cmd = f'grep -c "^>" {infile}'
    command = Command()
    command.command = cmd
    command.run(infile.resolve().parent)
    if command.stderr != '':
        raise RuntimeError(command.stderr, cmd)
    return int(command.stdout.rstrip())

def count_bases(infile: Path) -> int:
    """
    Counts the number of bases in a FASTA file.
    :param infile: Input FASTA file
    :return: Number of bases
    """
    total_size = 0
    for seq in SeqIO.parse(infile, 'fasta'):
        total_size += len(seq)
    return total_size

def gc(infile: Path) -> float:
    """
    Calculates the %GC-content across all sequences.
    The value is weighted by sequence length.
    :param infile: Input FASTA file
    :return: GC content (as percentage)
    """
    weights = []
    gcs = []
    for seq in SeqIO.parse(infile, 'fasta'):
        weights.append(len(seq))
        gcs.append(100 * gc_fraction(seq.seq))
    return float(np.average(gcs, weights=weights))

def is_indexed(fasta: Path) -> bool:
    """
    Checks if the input fasta file is indexed with the FAI index.
    :param fasta: Fasta file to check for associated index
    :return: True or False depending on the presence of the index file
    """
    try:
        next(fasta.parent.glob(f'{fasta.name}.fai'))
    except StopIteration:
        return False
    return True

def rename_sequences_regex(
        fasta_in: Path, fasta_out: Path, regex: str, repl: str, search_in_description: bool = False,
        sort_output: bool = True, description: Optional[str] = None) -> None:
    """
    Renames the sequences in the input FASTA file using a regex.
    :param fasta_in: Input FASTA file
    :param fasta_out: Output FASTA file
    :param regex: Regex
    :param repl: Replacement
    :param search_in_description: If True, the regex matches the full sequence description
    :param sort_output: If true, output records are sorted by sequence id
    :param description: Updated sequence description (optional)
    :return: None
    """
    with fasta_in.open() as handle:
        seqs = list(SeqIO.parse(handle, 'fasta'))
    for seq_record in seqs:
        seq_record.id = re.sub(regex, repl, seq_record.description if search_in_description else seq_record.id)
        if description is not None:
            seq_record.description = description
    if sort_output:
        seqs.sort(key=lambda x: x.id)
    with fasta_out.open('w') as handle:
        SeqIO.write(seqs, handle, 'fasta')

def rename_sequences_with_fasta_file(fasta_in: Path, fasta_ref: Path, fasta_out: Path) -> None:
    """
    Changes the names of the sequences in the input FASTA file to the sequences in the reference FASTA file (in the
    same order)
    :param fasta_in: Input FASTA file
    :param fasta_ref: Reference FASTA file
    :param fasta_out: Output FASTA file
    :return: None
    """
    # Parse reference FASTA file
    with fasta_ref.open() as handle:
        seq_ids_ref = [s.id for s in SeqIO.parse(handle, 'fasta')]

    # Parse input FASTA file
    with fasta_in.open() as handle:
        seqs_in = list(SeqIO.parse(handle, 'fasta'))

    # Update seq ids
    if len(seqs_in) != len(seq_ids_ref):
        raise ValueError(
            f'Nb. of seqs in input file ({len(seqs_in)}) does not match reference file ({len(seq_ids_ref)}')
    for seq_in, ref_id in zip(seqs_in, seq_ids_ref):
        seq_in.id = ref_id

    # Create the output file
    with fasta_out.open('w') as handle:
        SeqIO.write(seqs_in, handle, 'fasta')

def convert_fasta_to_fastq(fasta_input_file: Path, fastq_output_file: Path) -> None:
    """
    Converts a FASTA file to a FASTQ file with perfect Phred quality scores.
    :return: None
    """
    with fasta_input_file.open('r') as fasta_file, fastq_output_file.open('w') as fastq_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Create a fake quality score string of the same length as the sequence
            fake_quality = [40] * len(record.seq)

            # Create a SeqRecord with the same sequence and a fake quality string
            fake_record = SeqRecord(record.seq, id=record.id, description=record.description, letter_annotations={
                "phred_quality": fake_quality})
            # Write the SeqRecord in FASTQ format
            SeqIO.write(fake_record, fastq_file, 'fastq')

def has_duplicates(fasta: Path) -> bool:
    """
    Checks if there are duplicate sequence IDs in the FASTA file.
    :param fasta: FASTA file that has to be checked for duplicate sequence IDs
    :return: boolean indicating whether the FASTA file contains duplicate IDs or not
    """
    seq_ids = set()
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in seq_ids:
                return True
            seq_ids.add(record.id)
    return False
