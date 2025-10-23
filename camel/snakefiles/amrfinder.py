from pathlib import Path

DIR_AMRFINDER = Path('amrfinder')
SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Input
INPUT_FASTA = DIR_AMRFINDER / 'fasta.io'

# Report and summary
OUTPUT_TSV = DIR_AMRFINDER / 'tool' / 'tsv.io'
OUTPUT_INFORMS = DIR_AMRFINDER / 'tool' / 'informs.io'
OUTPUT_REPORT = DIR_AMRFINDER / 'report' / 'html.iob'
OUTPUT_REPORT_EMPTY = DIR_AMRFINDER / 'html-empty.iob'
OUTPUT_SUMMARY = DIR_AMRFINDER / 'summary_amrfinder.{ext}'
