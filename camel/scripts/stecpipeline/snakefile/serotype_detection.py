from pathlib import Path

SNAKEFILE_SEROTYPE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_VAL = 'serotype_detection/tool/val-sero.io'
OUTPUT_REPORT = 'serotype_detection/report/html.iob'
OUTPUT_SUMMARY = 'serotype_detection/summary/summary_out.{ext}'
