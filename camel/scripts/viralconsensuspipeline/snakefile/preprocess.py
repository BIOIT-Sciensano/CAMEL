from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_REPORT = 'preprocess/report/html.iob'
OUTPUT_SUMMARY = 'preprocess/report/summary_preprocess.{ext}'
OUTPUT_INFORMS = 'preprocess/report/informs.io'
OUTPUT_FASTQ = 'preprocess/downsample/fq_dict.io'

# AmpliGone
OUTPUT_AMPLIGONE_REPORT = 'preprocess/ampligone/report/html.iob'
OUTPUT_AMPLIGONE_REPORT_EMPTY = 'preprocess/ampligone/report/html-empty.iob'

# Amplicon clipping
OUTPUT_CLIPPING_REPORT = 'preprocess/ampliconclip/html.iob'
OUTPUT_CLIPPING_REPORT_EMPTY = 'preprocess/ampliconclip/html-empty.iob'
