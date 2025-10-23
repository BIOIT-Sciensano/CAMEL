from dataclasses import dataclass
from pathlib import Path

from camel.app.toolkits.blast.blasthitstatistics import BLASTN_OUTPUT_FORMAT
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.tools.blast.blastformatter import BlastFormatter
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.blast.blastx import Blastx
from camel.app.tools.kma.kma import KMA
from camel.app.tools.kma.kmatypinghitextractor import KMATypingHitExtractor
from camel.app.tools.pipelines.sequence_typing.alignmentextractor import AlignmentExtractor
from camel.app.tools.pipelines.sequence_typing.besthitselector import BestHitSelector


@dataclass(frozen=True, unsafe_hash=True)
class TypingResultHolder:
    """
    Class to store the result of a sequence typing job.
    """
    hit: ToolIOValue
    informs: dict[str, str]


def detect_hit_blast(
        dir_working: Path, fasta_in: Path, dir_scheme: Path, locus_metadata: dict, blastn_task: str = 'megablast',
        threads_per_job: int = 1) -> TypingResultHolder:
    """
    Performs hit detection with BLAST+.
    :param fasta_in: Input FASTA file
    :param dir_scheme: Base directory for the typing scheme
    :param locus_metadata: Metadata for the locus
    :param dir_working: Working directory
    :param blastn_task: blastn '-task' parameter value
    :param threads_per_job: Threads per BLAST job
    """
    # Create working directory
    dir_working = Path(dir_working)
    dir_working.mkdir(parents=True, exist_ok=True)

    # Initialize BLAST class
    if locus_metadata['type'] == 'DNA':
        blast = Blastn()
        blast.update_parameters(task=blastn_task)
    elif locus_metadata['type'] == 'peptide':
        blast = Blastx()
        blast.update_parameters(seg='no', comp_based_stats='0')
    else:
        raise ValueError(f"Invalid locus type: {locus_metadata['type']}")

    # Run BLAST+
    db_path = dir_scheme / locus_metadata['fasta_path']
    blast.add_input_files({'DB_BLAST': [ToolIOFile(db_path)], 'FASTA': [ToolIOFile(fasta_in)]})
    blast.update_parameters(threads=threads_per_job)
    blast.run(dir_working)

    # TSV generation
    formatter_tsv = BlastFormatter()
    formatter_tsv.update_parameters(output_format=BLASTN_OUTPUT_FORMAT)
    formatter_tsv.add_input_files({'ASN': blast.tool_outputs['ASN']})
    formatter_tsv.run(dir_working)

    # Best hit selection
    hit_selector = BestHitSelector()
    hit_selector.add_input_files({'TSV': formatter_tsv.tool_outputs['TSV']})
    hit_selector.add_input_informs({'locus': locus_metadata})
    hit_selector.run(dir_working)

    if not hit_selector.tool_outputs['VAL_Hit'][0].value.is_perfect_hit():
        # Text alignment generation
        formatter_text = BlastFormatter()
        formatter_text.update_parameters(output_format='0', num_alignments=1000)
        formatter_text.add_input_files({'ASN': blast.tool_outputs['ASN']})
        formatter_text.run(dir_working)

        # Alignment extraction
        extractor = AlignmentExtractor()
        extractor.add_input_files({
            'TXT': formatter_text.tool_outputs['TXT'], 'VAL_Hits': hit_selector.tool_outputs['VAL_Hit']})
        extractor.run(dir_working)

        # Add the alignment to the hit object
        if len(extractor.tool_outputs['TXT']) > 0:
            best_hit = hit_selector.tool_outputs['VAL_Hit'][0].value
            best_hit.alignment_path = extractor.tool_outputs['TXT'][0].path
    return TypingResultHolder(hit=hit_selector.tool_outputs['VAL_Hit'][0], informs=blast.informs)


# noinspection PyUnusedLocal
def detect_hit_kma(dir_working: Path, fastq_in: dict[str, list[ToolIOFile]], dir_scheme: Path, locus_metadata: dict,
                   read_type: str = 'illumina', threads_per_job: int = 1) -> TypingResultHolder:
    """
    Performs hit detection with KMA.
    :param fastq_in: Input FASTQ files
    :param dir_scheme: Base directory for the typing scheme
    :param locus_metadata: Metadata for the locus
    :param dir_working: Working directory
    :param threads_per_job: Threads per BLAST job
    :param read_type: Read type
    :return: Typing result
    """
    # Create working directory
    dir_working = Path(dir_working)
    dir_working.mkdir(parents=True, exist_ok=True)

    # Get the KMA database
    dir_locus = (dir_scheme / locus_metadata['fasta_path']).parent
    try:
        db_path = next((dir_locus / 'kma').glob('*.name'))
    except StopIteration:
        raise FileNotFoundError(f"KMA database for locus '{locus_metadata['name']}' ({dir_locus}) not found")

    # Launch KMA
    kma = KMA()
    kma.add_input_files(fastq_in)
    kma.add_input_files({'DB': [ToolIOValue(str(db_path.parent / db_path.stem))]})
    if read_type == 'nanopore':
        kma.update_parameters(bc_nano=None, basecalls='0.7')
    kma.run(dir_working)

    # Extract the best hit
    kma_extractor = KMATypingHitExtractor()
    kma_extractor.add_input_files({'TSV': kma.tool_outputs['TSV']})
    kma_extractor.add_input_informs({'locus': locus_metadata})
    kma_extractor.run(dir_working)
    return TypingResultHolder(hit=kma_extractor.tool_outputs['VAL_hit'][0], informs=kma.informs)


detection_by_method = {
    'blast': detect_hit_blast,
    'kma': detect_hit_kma,
}
