from pathlib import Path

from camel.app.snakemake import snakemakeutils
from camel.scripts.snpphylogeny.snakefile.trimming_all import TRIMMING_ALL

rule run_trimming_workflow:
    """
    This rule runs the trimming workflow on a sample.
    """
    input:
        FASTQ_PE = lambda wildcards: config['samples'][wildcards.sample]
    output:
        TRIMMING_OUTPUT = '{sample}/trimming/trimming_out.iob'
    threads: 4
    params:
        working_dir = lambda wildcards: wildcards.sample,
        adapter = config['read_trimming']['adapter']
    run:
        from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
        trimming = TrimmingIlluminaWrapper(Path(str(params.working_dir)).absolute() / 'trimming')
        trimming.run_workflow([Path(x) for x in input.FASTQ_PE], params.adapter, threads)
        snakemakeutils.dump_object(trimming.output, Path(output.TRIMMING_OUTPUT))

rule collect_trimming_output:
    """
    Combines the trimming output for all samples into a dictionary.
    """
    input:
        OUT_TRIM = expand(rules.run_trimming_workflow.output.TRIMMING_OUTPUT, sample=sorted(config['samples'].keys()))
    output:
        IO = TRIMMING_ALL
    params:
        samples = sorted(config['samples'].keys())
    run:
        output_data = {}
        for i, out in enumerate(input.OUT_TRIM):
            output_data[params.samples[i]] = snakemakeutils.load_object(Path(out))
        snakemakeutils.dump_object(output_data, Path(output.IO))
