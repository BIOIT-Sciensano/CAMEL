from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule assembly_spades_run:
    """
    De-novo assembly using SPAdes.
    """
    input:
        IO = 'fq_dict.io'
    output:
        FASTA_Contig = 'assembly/spades/fasta.io', # assembly_spades.OUTPUT_FASTA
        INFORMS = 'assembly/spades/informs.io' # assembly_spades.OUTPUT_INFORMS
    params:
        dir_ = 'assembly/spades',
        spades_options = config.get('assembly', {}).get('spades', {})
    threads: 8
    priority: 1
    run:
        from camel.app.tools.spades.spades import SPAdes
        from camel.app.core.snakemake import snakepipelineutils
        spades = SPAdes()

        # Reformat FASTQ dictionary
        fq_dict = snakepipelineutils.extract_fq_input(
            Path(input.IO),
            key_pe='FASTQ_PE_1',
            keys_se=['FASTQ_SE_1', 'FASTQ_SE_2'],
            key_se='FASTQ_SE_1',
            drop_empty=True,
            read_type='PE')
        spades.add_input_files(fq_dict)
        step = Step(rule_name=str(rule), tool=spades, dir_=Path(params.dir_))
        spades.update_parameters(**params.spades_options)
        spades.update_parameters(isolate=True, careful=False, threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(spades, output)
