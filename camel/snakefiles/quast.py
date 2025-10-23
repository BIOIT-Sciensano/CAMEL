from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'quast/report/html.iob'
OUTPUT_SUMMARY = 'quast/summary/summary_quast.{ext}'
OUTPUT_INFORMS = 'quast/output/informs.io'
OUTPUT_INFORMS_BUSCO = 'quast/busco/informs.io'
