from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Input
INPUT_FASTA = 'btyper/input/fasta.io'

# Report and summary
OUTPUT_VAL = 'btyper/tool/val-btyper.io'
OUTPUT_INFORMS = 'btyper/tool/informs.io'
OUTPUT_REPORT = 'btyper/report/html.iob'
OUTPUT_REPORT_EMPTY = 'btyper/report/html-empty.iob'
OUTPUT_SUMMARY = 'btyper/summary/summary_out.{ext}'
