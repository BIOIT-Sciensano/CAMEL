from pathlib import Path

from camelcore.app.command import Command

from camel.app.core.errors import SnakemakeExecutionError
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.broadwgs import references
from camel.scripts.broadwgs.snakefile import alignment, bam_to_cram


rule samtools_convert_to_cram:
    """
    Convert the aligned BAM to CRAM for storage.
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
    output:
        CRAM = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert",
        output_file = f'{config["sample"]}.cram'
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bam_to_cram = SamtoolsView()
        snakemakeutils.add_io_input(bam_to_cram,"BAM", Path(input.BAM))
        snakemakeutils.add_io_input(bam_to_cram,"FASTA_REF", Path(input.FASTA_REF))
        step = Step(rule_name=str(rule), tool=bam_to_cram, dir_=params.working_dir)
        bam_to_cram.update_parameters(
            output_filename = params.output_file,
            output_format = 'CRAM',
            threads = threads,
            **config['rule_params']['bam_to_cram'][rule]
        )
        step.run()
        snakemakeutils.dump_io_output(bam_to_cram,'CRAM', Path(output.CRAM))

rule checksum_cram:
    """
    Calculate the checksum of the CRAM file.
    """
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
    output:
        CRAM_checksum = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM_checksum
    threads: config["params_smk"]["threads_cram"]
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert"
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        cram_file = snakemakeutils.load_object(Path(input.CRAM))[0]

        cmd_checksum = Command(f"md5sum {cram_file} > {output.CRAM_checksum}")
        cmd_checksum.run(params.working_dir)
        if cmd_checksum.exit_code != 0:
            raise SnakemakeExecutionError(cmd_checksum.stdout, cmd_checksum.stderr)

rule samtools_index_cram:
    """
    Index the CRAM file.
    """
    input:
        CRAM = rules.samtools_convert_to_cram.output.CRAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE
    output:
        CRAI = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAI,
    params:
        working_dir = Path(config['working_dir']) / "bamtocram" / "convert",
        output_file = f'{config["sample"]}.cram'
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram

        Path(params.working_dir).mkdir(exist_ok=True)

        cram_index = SamtoolsIndexCram()
        snakemakeutils.add_io_inputs(cram_index, input)
        step = Step(rule_name=str(rule), tool=cram_index, dir_=params.working_dir)
        cram_index.update_parameters(output_filename = params.working_dir / params.output_file)
        step.run()
        snakemakeutils.dump_io_output(cram_index,'CRAI', Path(output.CRAI))
