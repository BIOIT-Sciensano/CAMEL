from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import variant_calling_clair3


rule clair3_variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        FASTA = variant_calling_clair3.OUTPUT_FASTA,
        INFORMS = 'variant_calling/reference/informs.io'
    params:
        reference = config['variant_calling'].get('reference')
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        snakemakeutils.dump_object([ToolIOFile(Path(params.reference['path']))], Path(output.FASTA))
        snakemakeutils.dump_object(params.reference, Path(output.INFORMS))

rule clair3_variant_calling_alignment_sorting:
    """
    Sorts the alignment.
    """
    input:
        BAM = variant_calling_clair3.INPUT_BAM
    output:
        BAM = variant_calling_clair3.OUTPUT_BAM
    params:
        running_dir = 'variant_calling/alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort()
        step = Step(rule_name=str(rule), tool=samtools_sort, dir_=Path(str(params.running_dir)))
        snakemakeutils.add_io_inputs(samtools_sort, input)
        step.run()
        snakemakeutils.dump_io_outputs(samtools_sort, output)

rule clair3_variant_calling_index_bam:
    """
    Index the bam file.
    """
    input:
        BAM = rules.clair3_variant_calling_alignment_sorting.output.BAM
    output:
        BAM = variant_calling_clair3.OUTPUT_BAM_INDEX
    params:
        running_dir = 'variant_calling/alignment_sorting'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex()
        step = Step(rule_name=str(rule), tool=samtools_index, dir_=Path(str(params.running_dir)))
        snakemakeutils.add_io_inputs(samtools_index, input)
        step.run()
        snakemakeutils.dump_io_outputs(samtools_index, output)

rule clair3_variant_calling_index_fasta:
    """
    Index the fasta file.
    """
    input:
        FASTA = rules.clair3_variant_calling_prep_reference.output.FASTA
    output:
        FASTA = variant_calling_clair3.OUTPUT_FASTA_INDEX
    params:
        running_dir = 'variant_calling/reference'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools_fastaindex = SamtoolsFastaIndex()
        step = Step(rule_name=str(rule), tool=samtools_fastaindex, dir_=Path(str(params.running_dir)))
        snakemakeutils.add_io_inputs(samtools_fastaindex, input)
        step.run()
        snakemakeutils.dump_io_outputs(samtools_fastaindex, output)

rule clair3_variant_calling:
    """
    This step actually calls variants with clair3
    """
    input:
        FASTA = rules.clair3_variant_calling_index_fasta.output.FASTA,
        BAM = rules.clair3_variant_calling_alignment_sorting.output.BAM,
        BAM_INDEX = rules.clair3_variant_calling_index_bam.output.BAM
    output:
        VCF = 'variant_calling/calling/vcf_gz.io',
        INFORMS = 'variant_calling/calling/informs.io'
    params:
        running_dir = 'variant_calling/calling',
        long_indel = config['variant_calling'].get('clair3', {}).get('long_indel', False),
        haploid_precise = config['variant_calling'].get('clair3', {}).get('haploid_precise', False),
        include_ctgs = config['variant_calling'].get('clair3', {}).get('include_ctgs', False),
        no_phasing = config['variant_calling'].get('clair3', {}).get('no_phasing', False),
        model_path = config['variant_calling'].get('clair3', {}).get('model_path'),
        platform = config['variant_calling'].get('clair3', {}).get('platform')
    threads: 8
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        clair3 = Clair3()
        snakemakeutils.add_io_inputs(clair3,input)
        step = Step(rule_name=str(rule), tool=clair3, dir_=Path(str(params.running_dir)))

        # Update the parameters
        if params.model_path is None:
            raise ValueError('"model_path" is required')
        if params.platform is None:
            raise ValueError('"platform" is required')
        clair3.update_parameters(
            model_path=params.model_path,
            long_indel=params.long_indel,
            include_ctgs=params.include_ctgs,
            no_phasing=params.no_phasing,
            platform=params.platform,
            output_path=str(params.running_dir),
            threads=threads
        )

        # Run the tool
        step.run()
        snakemakeutils.dump_io_outputs(clair3, output)

rule clair3_variant_calling_normalize_indels:
    """
    Normalizes indels.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling.output.VCF,
        FASTA = rules.clair3_variant_calling_prep_reference.output.FASTA
    output:
        VCF_GZ = 'variant_calling/norm/vcf_gz.io'
    params:
        running_dir = 'variant_calling/norm'
    run:
        from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm
        bcftools_norm = BcftoolsNorm()
        snakemakeutils.add_io_inputs(bcftools_norm, input)
        step = Step(rule_name=str(rule), tool=bcftools_norm, dir_=Path(str(params.running_dir)))
        bcftools_norm.update_parameters(output_type='z')
        step.run()
        snakemakeutils.dump_io_outputs(bcftools_norm, output)

rule clair3_variant_calling_index_vcf_gz:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling_normalize_indels.output.VCF_GZ
    output:
        VCF_GZ = variant_calling_clair3.OUTPUT_UNFILTERED_VCF_GZ
    params:
        running_dir = 'variant_calling/norm'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        indexer = BcftoolsIndex()
        snakemakeutils.add_io_inputs(indexer, input)
        step = Step(rule_name=str(rule), tool=indexer, dir_=Path(str(params.running_dir)))
        step.run()
        snakemakeutils.dump_io_outputs(indexer, output)

rule clair3_variant_calling_unzip_vcf:
    """
    Unzips the VCF file.
    """
    input:
        VCF_GZ = rules.clair3_variant_calling_index_vcf_gz.output.VCF_GZ
    output:
        VCF = variant_calling_clair3.OUTPUT_UNFILTERED_VCF
    params:
        running_dir = 'variant_calling/unzip_vcf',
        sample_name = config['input']['sample_name']
    run:
        from camel.app.core.utils import fileutils
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView()
        snakemakeutils.add_io_inputs(bcftools_view, input)
        step = Step(rule_name=str(rule), tool=bcftools_view, dir_=Path(str(params.running_dir)))
        output_filename = f'variants-{fileutils.make_valid(params.sample_name)}.vcf'
        bcftools_view.update_parameters(output_type='v', output_filename=output_filename)
        step.run()
        snakemakeutils.dump_io_outputs(bcftools_view, output)
