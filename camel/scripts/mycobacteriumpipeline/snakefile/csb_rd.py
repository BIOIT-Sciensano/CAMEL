from pathlib import Path

SNAKEFILE_CSB_RD = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_CSB_RD_REPORT = 'csb_rd/report/html.iob'
OUTPUT_CSB_RD_REPORT_EMPTY = 'csb_rd/report/html-empty.iob'
OUTPUT_CSB_RD_SUMMARY = 'csb_rd/summary.{ext}'
