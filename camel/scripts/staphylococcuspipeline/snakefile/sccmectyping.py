from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'sccmec_typing/report/html.iob'
OUTPUT_REPORT_EMPTY = 'sccmec_typing/report/html-empty.iob'
OUTPUT_SUMMARY = 'sccmec_typing/summary_sccmec.{ext}'
