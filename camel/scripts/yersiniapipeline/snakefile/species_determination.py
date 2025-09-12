from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'species_determination/report/html.iob'
OUTPUT_REPORT_EMPTY = 'species_determination/report/html-empty.iob'
OUTPUT_INFORMS = 'species_determination/tool/informs.io'
OUTPUT_SUMMARY = 'species_determination/summary/summary.{ext}'
