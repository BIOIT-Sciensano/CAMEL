from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_FASTA = 'assembly/flye/fasta.io'
OUTPUT_INFORMS = 'assembly/flye/informs.io'
