from pathlib import Path

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

INPUT_ASSEMBLY_FASTA = 'polish/short_reads/{assembly_type}/input/fasta.io'
OUTPUT_POLISHING_FASTA = 'polish/short_reads/{assembly_type}/polca/fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLYPOLISH = 'polish/short_reads/{assembly_type}/polypolish/fasta-index.io'
OUTPUT_POLYPOLISH_INFORMS = 'polish/short_reads/{assembly_type}/polypolish/informs.io'
OUTPUT_POLYPOLISH_FASTA = 'polish/short_reads/{assembly_type}/polypolish/fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLCA = 'polish/short_reads/{assembly_type}/polca/fasta-index.io'
OUTPUT_POLCA_INFORMS = 'polish/short_reads/{assembly_type}/polca/informs.io'
