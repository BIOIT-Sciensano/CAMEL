from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'amr/report/html.iob'
OUTPUT_REPORT_CDS = 'amr/cds/html.iob'
OUTPUT_REPORT_EMPTY = 'amr/report/html-empty.iob'
OUTPUT_SUMMARY = 'amr/summary/summary_out.{ext}'
OUTPUT_INFORMS_CSQ = 'amr/csq/informs.io'
