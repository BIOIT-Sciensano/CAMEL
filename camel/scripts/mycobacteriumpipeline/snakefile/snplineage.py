from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'snp_lineage/report/html.iob'
OUTPUT_REPORT_EMPTY = 'snp_lineage/report/html-empty.io'
OUTPUT_SUMMARY = 'snp_lineage/summary/summary_out.{ext}'
OUTPUT_INFORMS = 'snp_lineage/tool/informs.iob'
