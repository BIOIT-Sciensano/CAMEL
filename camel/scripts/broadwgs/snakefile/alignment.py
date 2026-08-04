from pathlib import Path

SNAKEFILE_ALIGNMENT = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_alignment = Path('alignment')
_dir_qc = Path('qc')
OUTPUT_INTERMEDIATE_BAM = _dir_alignment/ "add_readgroups" / "{input_basename}.aligned_rgadded.bam.io"
OUTPUT_ALIGNMENT_BAM = _dir_alignment / "gather_bqsr_sorted_bam" / "bqsr_gathered_sorted.bam.io"
OUTPUT_MARK_DUPLICATES_METRICS = _dir_qc / 'mark_duplicates' / "duplicate_metrics.txt.io"
