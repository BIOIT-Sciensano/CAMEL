from pathlib import Path

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

OUTPUT_REPORT = 'ref_selection/report/html.iob'
OUTPUT_REPORT_EMPTY = 'ref_selection/report/html-empty.iob'
OUTPUT_FASTA = 'ref_selection/create_fasta/fasta.io'
OUTPUT_SUMMARY = 'ref_selection/summary/summary.{ext}'
