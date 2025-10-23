from pathlib import Path

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

INPUT_ASSEMBLY_FASTA = 'polish/short_reads/{assembly_type}/input/fasta.io'
OUTPUT_POLISHING_FASTA = 'polish/short_reads/{assembly_type}/pypolca/fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLYPOLISH = 'polish/short_reads/{assembly_type}/polypolish/fasta-index.io'
OUTPUT_POLYPOLISH_INFORMS = 'polish/short_reads/{assembly_type}/polypolish/informs.io'
OUTPUT_POLYPOLISH_FASTA = 'polish/short_reads/{assembly_type}/polypolish/fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_PYPOLCA = 'polish/short_reads/{assembly_type}/pypolca/fasta-index.io'
OUTPUT_PYPOLCA_INFORMS = 'polish/short_reads/{assembly_type}/pypolca/informs.io'
