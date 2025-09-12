from pathlib import Path

from camel.resources.snakefile import assembly_flye

SNAKEFILE_QC = f'{Path(__file__).parent / Path(__file__).stem}.smk'

consensus_by_tool = {
    'Flye':  assembly_flye.OUTPUT_FASTA,
    'Medaka': 'polish/long_reads/flye/sequence/fasta.io',
    'Polypolish': 'polish/short_reads/flye/polypolish/fasta.io',
    'POLCA': 'polish/short_reads/flye/polca/fasta.io',
    'Unicycler': 'unicycler/fasta.io',
    'Medaka-Unicycler': 'polish/long_reads/unicycler/sequence/fasta.io',
    'Polypolish-Unicycler': 'polish/short_reads/unicycler/polypolish/fasta.io',
    'POLCA-Unicycler': 'polish/short_reads/unicycler/polca/fasta.io'
}

ALE_KEYS = ['depth', 'kmer', 'insert', 'place']
