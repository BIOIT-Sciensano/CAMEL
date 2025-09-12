from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Report and summary
OUTPUT_INFORMS = 'kleborate/tool/informs.io'
OUTPUT_REPORT = 'kleborate/report/html.iob'
OUTPUT_REPORT_EMPTY = 'kleborate/report/html-empty.iob'
OUTPUT_SUMMARY = 'kleborate/summary/summary_kleborate.{ext}'
