from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Report and summary
OUTPUT_INFORMS = 'bacmet/blastp/informs.io'
OUTPUT_REPORT = 'bacmet/report/html.iob'
OUTPUT_REPORT_EMPTY = 'bacmet/report/html-empty.iob'
OUTPUT_SUMMARY = 'bacmet/summary_bacmet.tsv'

# Prodigal
OUTPUT_PRODIGAL_INFORMS = 'bacmet/prodigal/tool/informs.io'
OUTPUT_PRODIGAL_REPORT = 'bacmet/prodigal/report/html.iob'
OUTPUT_PRODIGAL_REPORT_EMPTY = 'bacmet/prodigal/report/html-empty.iob'
