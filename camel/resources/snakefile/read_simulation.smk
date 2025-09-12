from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import assembly

rule read_simulation_fastq_from_fasta:
    """
    Simulates illumina reads from the input FASTA file.
    This rule is executed when the input type is FASTA and no VCF input file has been provided.
    """
    input:
        FASTA = assembly.get_fasta_raw(config)
    output:
        FASTQ_PE = 'read_simulation/art/fastq.io', # read_simulation.OUTPUT_FASTQ
        INFORMS = 'read_simulation/art/informs.io' # read_simulation.OUTPUT_INFORMS
    params:
        dir_ = 'read_simulation/art'
    run:
        from camel.app.tools.art.art import ART
        art = ART()
        step = Step(rule_name=str(rule), tool=art, dir_=Path(str(params.dir_)))
        snakemakeutils.add_pickle_inputs(art, input)
        step.run()
        snakemakeutils.dump_tool_outputs(art, output)
