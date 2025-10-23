from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'gmats/report/html.iob'
OUTPUT_REPORT_EMPTY = 'gmats/report/html-empty.iob'
OUTPUT_SUMMARY = 'gmats/summary/summary_gmats.{ext}'
