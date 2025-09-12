from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'antivirals/report/html.iob'
OUTPUT_REPORT_EMPTY = 'antivirals/report/html-empty.iob'
OUTPUT_SUMMARY = 'antivirals/summary/summary.{ext}'
