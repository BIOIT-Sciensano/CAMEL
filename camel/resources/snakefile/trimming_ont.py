from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

INPUT_ONT_FASTQ = 'trimming_ont/input/fastq.io'

# Report and summary
OUTPUT_REPORT = 'trimming_ont/report/html.iob'
OUTPUT_SUMMARY = 'trimming_ont/summary/summary_out.{ext}'

# Trimming
OUTPUT_READS = 'trimming_ont/seqkit/fastq.io'
OUTPUT_INFORMS = 'trimming_ont/seqkit/informs.io'
OUTPUT_DICT = 'trimming_ont/fastq_all.io'

# NanoPlot
OUTPUT_NANOPLOT_TXT_PRE = 'trimming_ont/nanoplot-pre/txt.io'
OUTPUT_NANOPLOT_HTML_PRE = 'trimming_ont/nanoplot-pre/html.io'
OUTPUT_NANOPLOT_INFORMS_PRE = 'trimming_ont/nanoplot-pre/informs.io'
OUTPUT_NANOPLOT_TXT_POST = 'trimming_ont/nanoplot-post/txt.io'
OUTPUT_NANOPLOT_HTML_POST = 'trimming_ont/nanoplot-post/html.io'
OUTPUT_NANOPLOT_INFORMS_POST = 'trimming_ont/nanoplot-post/informs.io'
