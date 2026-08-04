from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_INFORMS = 'spa_typing/blastn/informs.io'
OUTPUT_REPORT = 'spa_typing/report/html.iob'
OUTPUT_REPORT_EMPTY = 'spa_typing/report/html-empty.iob'
OUTPUT_SUMMARY = 'spa_typing/summary_spatyping.{ext}'
