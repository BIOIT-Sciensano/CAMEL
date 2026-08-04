from pathlib import Path

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

INPUT_ASSEMBLY_FASTA = 'polish/long_reads/{assembly_type}/input/fasta.io'
OUTPUT_FASTA = 'polish/long_reads/{assembly_type}/sequence/fasta.io'
OUTPUT_POLISH_MEDAKA_INFORMS = 'polish/long_reads/{assembly_type}/sequence/commands-sequence.io'
OUTPUT_ASSEMBLY_REPORT_EMPTY = 'polish/long_reads/{assembly_type}/report/html-empty.iob'
