from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'snpit/html.iob'
OUTPUT_REPORT_EMPTY = 'snpit/html-empty.iob'
OUTPUT_SUMMARY = 'snpit/summary_out.{ext}'
OUTPUT_INFORMS = 'snpit/informs.io'
