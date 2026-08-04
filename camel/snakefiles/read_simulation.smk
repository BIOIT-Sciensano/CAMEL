from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly

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
    run:
        from camel.app.tools.art.art import ART
        art = ART()
        step = Step(rule_name=str(rule), tool=art, dir_=snakemakeutils.get_rule_dir(output))
        snakemakeutils.add_io_inputs(art, input)
        step.run()
        snakemakeutils.dump_io_outputs(art, output)
