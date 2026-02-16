from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_UNFILTERED_VCF = 'variant_calling/vcf.io'
OUTPUT_REPORT_LOFREQ = 'variant_calling/report/html-lofreq.iob'
OUTPUT_LOFREQ_INFORMS = 'variant_calling/informs.io'
