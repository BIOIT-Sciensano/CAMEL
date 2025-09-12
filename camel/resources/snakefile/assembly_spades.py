from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_FASTA = 'assembly/spades/fasta.io'
OUTPUT_INFORMS = 'assembly/spades/informs.io'
