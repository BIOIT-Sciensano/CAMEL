from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'spoligotyping/report/html.iob'
OUTPUT_REPORT_EMPTY = 'spoligotyping/report/html-empty.iob'
OUTPUT_SUMMARY = 'spoligotyping/summary/summary_out.{ext}'
OUTPUT_INFORMS = 'spoligotyping/spotyping/informs.io'
