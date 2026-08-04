from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'multi_allelic/report/html.iob'
OUTPUT_SUMMARY = 'multi_allelic/report/summary.{ext}'
