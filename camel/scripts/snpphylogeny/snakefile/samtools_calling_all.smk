from pathlib import Path

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.snpphylogeny.snakefile.samtools_calling_all import OUTPUT_CALLING_ALL

rule run_variant_calling_workflow:
    """
    This rule runs the variant calling workflow on a sample.
    """
    output:
        IO = '{sample}/variant_calling_out.iob'
    threads: 4
    params:
        working_dir = lambda wildcards: wildcards.sample,
        ref_info = config['reference_info'],
        calling_options = config['options'],
        sample_name = lambda wildcards: wildcards.sample,
        sample_config = lambda wildcards: config['samples'][wildcards.sample]
    run:
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput
        from camel.app.wrappers.variantcallingwrapper import VariantCallingWrapper
        fastq_input = FastqInput(
            'illumina',
            pe=[ToolIOFile(Path(x)) for x in params.sample_config['PE']],
            se_fwd=[ToolIOFile(Path(params.sample_config['SE_FWD']))] if 'SE_FWD' in params.sample_config else None,
            se_rev=[ToolIOFile(Path(params.sample_config['SE_REV']))] if 'SE_REV' in params.sample_config else None
        )
        wrapper = VariantCallingWrapper(Path(str(params.working_dir)).absolute() / 'variant_calling')
        wrapper.run(
            params.ref_info, str(params.sample_name), fastq_input, 'illumina', params.calling_options, int(str(threads)))
        snakemakeutils.dump_object(wrapper.output, Path(output.IO))

rule collect_variant_calling_output:
    """
    Combines the variant calling output for all samples into a dictionary.
    """
    input:
        VC_OUT = expand(rules.run_variant_calling_workflow.output.IO, sample=sorted(list(config['samples'].keys())))
    output:
        IO = OUTPUT_CALLING_ALL
    params:
        samples = sorted(config['samples'].keys())
    run:
        output_data = {}
        for i, out in enumerate(input.VC_OUT):
            output_data[params.samples[i]] = snakemakeutils.load_object(Path(out))
        snakemakeutils.dump_object(output_data, Path(output.IO))
