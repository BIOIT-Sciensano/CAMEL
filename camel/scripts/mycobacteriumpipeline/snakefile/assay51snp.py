from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = '51snp/report/html.iob'
OUTPUT_REPORT_EMPTY = '51snp/report/html-empty.iob'
OUTPUT_SUMMARY = '51snp/summary_out.{ext}'
