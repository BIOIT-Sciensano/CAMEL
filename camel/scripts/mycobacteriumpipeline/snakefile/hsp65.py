from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'hsp65/report/html.iob'
OUTPUT_REPORT_EMPTY = 'hsp65/report/html-empty.iob'
OUTPUT_SUMMARY = 'hsp65/summary/summary.{ext}'
