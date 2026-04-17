from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule assembly_flye_run:
    """
    De-novo assembly using Flye.
    """
    input:
        IO = 'fq_dict.io'
    output:
        FASTA = 'assembly/flye/fasta.io', # assembly_flye.OUTPUT_FASTA
        INFORMS = 'assembly/flye/informs.io' # assembly_flye.OUTPUT_INFORMS
    params:
        dir_ = 'assembly/flye',
        flye_options = config.get('assembly', {}).get('flye', {})
    threads: 8
    priority: 1
    run:
        from camel.app.tools.flye.flye import Flye
        from camel.app.core.snakemake import snakepipelineutils
        flye = Flye()

        # Reformat FASTQ dictionary
        fq_dict = snakepipelineutils.extract_fq_input(Path(input.IO), key_se='FASTQ', read_type='SE')
        flye.add_input_files(fq_dict)
        step = Step(rule_name=str(rule), tool=flye, dir_=Path(params.dir_))
        flye.update_parameters(**params.flye_options)
        flye.update_parameters(threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(flye, output)
