from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# FASTA input
INPUT_FASTA = 'ani/input/fasta.io'

# Tool outputs
OUTPUT_VAL = 'ani/tool/val-ani.io'
OUTPUT_INFORMS = 'ani/tool/informs.io'

# Report and summary
OUTPUT_REPORT = 'ani/report/html.iob'
OUTPUT_REPORT_EMPTY = 'ani/report/html-empty.iob'
OUTPUT_SUMMARY = 'ani/summary/summary_out.{ext}'
