from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_FASTQ = 'read_simulation/art/fastq.io'
OUTPUT_INFORMS = 'read_simulation/art/informs.io'
