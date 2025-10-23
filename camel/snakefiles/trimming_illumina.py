from pathlib import Path

_dir = Path('trimming_illumina')
SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Input
INPUT_FASTQ = 'trimming_illumina/input/fastq.io'

# Report and summary
OUTPUT_REPORT = 'trimming_illumina/report/html.iob'
OUTPUT_SUMMARY = 'trimming_illumina/summary/summary_trim.{ext}'

# Trimming
OUTPUT_READS_PE = 'trimming_illumina/reads/fastq-pe.io'
OUTPUT_READS_SE_FWD = 'trimming_illumina/reads/fastq-se-fwd.io'
OUTPUT_READS_SE_REV = 'trimming_illumina/reads/fastq-se-rev.io'
OUTPUT_DICT = 'trimming_illumina/reads/fq_dict.io'
OUTPUT_INFORMS = 'trimming_illumina/report/informs.io'

# FastQC
OUTPUT_FASTQC_TXT_PRE = 'trimming_illumina/fastqc-pre/txt.io'
OUTPUT_FASTQC_HTML_PRE = 'trimming_illumina/fastqc-pre/html.io'
OUTPUT_FASTQC_TXT_POST = 'trimming_illumina/fastqc-post/txt.io'
OUTPUT_FASTQC_HTML_POST = 'trimming_illumina/fastqc-post/html.io'


def select_fastq_output(config: dict) -> Path:
    """
    Selects the PE FASTQ output based on the Snakemake configuration.
    :param config: Config data
    :return: Path to PE FASTQ output
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'trimmomatic' / 'fastq-pe.io'
    else:
        return Path('trimming_illumina') / 'fastp' / 'fastq-pe.io'


def select_informs(config: dict) -> Path:
    """
    Selects the trimming informs based on the Snakemake configuration.
    :param config: Config data
    :return: Path to the trimming informs
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'trimmomatic' / 'informs.io'
    else:
        return Path('trimming_illumina') / 'fastp' / 'informs.io'


def select_report_output(config: dict) -> Path:
    """
    Selects the trimming report based on the Snakemake configuration.
    :param config: Config data
    :return: Path to the output report
    """
    method = config['read_trimming'].get('method', 'trimmomatic')
    if method == 'trimmomatic':
        return Path('trimming_illumina') / 'report' / 'trimmomatic' / 'html.iob'
    else:
        return Path('trimming_illumina') / 'report' / 'fastp' / 'html.iob'
