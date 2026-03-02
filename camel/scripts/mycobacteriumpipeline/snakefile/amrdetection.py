from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = 'amr/report/html.iob'
OUTPUT_REPORT_CDS = 'amr/cds/html.iob'
OUTPUT_REPORT_EMPTY = 'amr/report/html-empty.iob'
OUTPUT_SUMMARY = 'amr/summary/summary_out.{ext}'
OUTPUT_INFORMS_CSQ = 'amr/csq/{caller}/informs.io'
OUTPUT_INFORMS_LOFREQ = 'amr/lofreq/vcf/informs.io'

OUTPUT_INFORMS = [OUTPUT_INFORMS_CSQ.format(caller='bcftools'), OUTPUT_INFORMS_LOFREQ]
