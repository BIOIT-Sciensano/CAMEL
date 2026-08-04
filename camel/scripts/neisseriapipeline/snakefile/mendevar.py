from pathlib import Path


SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'mendevar/report/html.iob'
OUTPUT_REPORT_EMPTY = 'mendevar/report/html-empty.iob'
OUTPUT_SUMMARY = 'mendevar/summary/summary_mendevar.{ext}'
