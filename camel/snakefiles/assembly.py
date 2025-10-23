from pathlib import Path
from typing import Any

from camel.snakefiles import assembly_spades, assembly_flye, polish_assembly_short, polish_assembly_long, \
    human_read_scrubbing

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_FASTA = 'assembly/filtering/fasta.io'
OUTPUT_INFORMS_FILTERING = 'assembly/filtering/informs.io'
OUTPUT_REPORT = 'assembly/report/html.iob'
OUTPUT_SUMMARY = 'assembly/summary/summary_out.tsv'
OUTPUT_INFORMS_MAPPING = 'assembly/{mapper}/informs.io'
OUTPUT_INFORMS_DEPTH = 'assembly/samtools_depth/{mapper}/informs.io'
OUTPUT_INFORMS_MAPPING_RATE = 'assembly/samtools_flagstat/{mapper}/informs.io'


def get_fasta_raw(config: dict[str, Any]) -> str:
    """
    Returns the assembly FASTA output IO object path (before filtering).
    """
    if (config['input_type'] in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' not in config['analyses']):
        return human_read_scrubbing.INPUT_FASTA.format(input_format='fasta')
    if (config['input_type'] in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        return human_read_scrubbing.OUTPUT_FASTA.format(input_format='fasta')
    if config['input_type'] == 'illumina':
        return assembly_spades.OUTPUT_FASTA
    if config['input_type'] == 'ont':
        return assembly_flye.OUTPUT_FASTA
    if config['input_type'] == 'hybrid':
        return polish_assembly_short.OUTPUT_POLISHING_FASTA.format(assembly_type='flye')
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_command_informs(config: dict[str, Any]) -> list[str]:
    """
    Returns the assembly informs output IO object paths.
    :return: Assembly informs path
    """
    if config['input_type'] in ('fasta', 'fasta_with_vcf'):
        return []
    if config['input_type'] == 'illumina':
        return [assembly_spades.OUTPUT_INFORMS]
    if config['input_type'] == 'ont':
        return [assembly_flye.OUTPUT_INFORMS]
    if config['input_type'] == 'hybrid':
        return [
            assembly_flye.OUTPUT_INFORMS,
            polish_assembly_long.OUTPUT_POLISH_MEDAKA_INFORMS.format(assembly_type='flye'),
            polish_assembly_short.OUTPUT_POLYPOLISH_INFORMS.format(assembly_type='flye'),
            polish_assembly_short.OUTPUT_PYPOLCA_INFORMS.format(assembly_type='flye')
        ]
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_mapping_inform(read_key: str) -> Path:
    """
    Returns the mapping informs.
    :param read_key: Read key
    :return: Path to mapping informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_INFORMS_MAPPING).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_INFORMS_MAPPING).format(mapper='minimap2'))
    raise ValueError(f'Invalid read key: {read_key}')


def get_depth_inform(read_key: str) -> Path:
    """
    Returns the depth informs.
    :param read_key: Read key
    :return: Path to depth informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_INFORMS_DEPTH).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_INFORMS_DEPTH).format(mapper='minimap2'))
    raise ValueError(f'Invalid read key: {read_key}')


def get_mapping_rate_inform(read_key: str) -> Path:
    """
    Returns the mapping rate informs.
    :param read_key: Read key
    :return: Path to depth informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_INFORMS_MAPPING_RATE).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_INFORMS_MAPPING_RATE).format(mapper='minimap2'))
    raise ValueError(f'Invalid read key: {read_key}')


def get_qc_informs(input_type: str, mode: str = 'assembly') -> list[Path]:
    """
    Returns the QC informs based on the input type.
    :param input_type: Input type
    :param mode: Reference mode (assembly or ref)
    :return: List of paths to QC informs
    """
    informs = []
    if input_type in ('hybrid', 'illumina'):
        informs.append(get_mapping_inform('fastq_pe'))
        if mode == 'assembly':
            informs.append(get_depth_inform('fastq_pe'))
    if input_type in ('hybrid', 'ont'):
        informs.append(get_mapping_inform('fastq_se'))
        if mode == 'assembly':
            informs.append(get_depth_inform('fastq_se'))
    return informs
