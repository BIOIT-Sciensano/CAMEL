from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import polish_assembly_short


rule polishing_samtools_index_polypolish:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA = lambda wildcards: polish_assembly_short.INPUT_ASSEMBLY_FASTA.format(assembly_type=wildcards.assembly_type)
    output:
        FASTA = 'polish/short_reads/{assembly_type}/polypolish/fasta-index.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/polypolish'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex()
        snakemakeutils.add_pickle_inputs(samtools, input)
        step = Step(rule_name=str(rule), tool=samtools, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools, output)

rule polishing_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.polishing_samtools_index_polypolish.output.FASTA
    output:
        INDEX_GENOME_PREFIX = 'polish/short_reads/{assembly_type}/polypolish/genome_prefix.iob'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/polypolish'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex()
        snakemakeutils.add_pickle_inputs(bwa, input)
        step = Step(rule_name=str(rule), tool=bwa, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bwa, output)

rule polishing_read_mapping_1:
    """
    Maps the forward reads against the assembly.
    """
    input:
        FQ_dict = 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX_SAMTOOLS = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = 'polish/short_reads/{assembly_type}/read_mapping/forward/bwa_readmap.io',
        INFORMS = 'polish/short_reads/{assembly_type}/read_mapping/forward/commands.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/read_mapping/forward'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap()
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[0]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        snakemakeutils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(rule_name=str(rule), tool=bwa_map, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bwa_map, output)

rule polishing_read_mapping_2:
    """
    Maps the reverse reads against the assembly.
    """
    input:
        FQ_dict = 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA,
        INDEX_GENOME_PREFIX = rules.polishing_bwa_index.output.INDEX_GENOME_PREFIX
    output:
        SAM = 'polish/short_reads/{assembly_type}/read_mapping/reverse/bwa_readmap.io',
        INFORMS = 'polish/short_reads/{assembly_type}/read_mapping/reverse/commands.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/read_mapping/reverse'
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap()
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_SE': [fq_in.pe[1]]})
        bwa_map.update_parameters(threads=threads, all_alns=True)
        snakemakeutils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))
        step = Step(rule_name=str(rule), tool=bwa_map, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bwa_map, output)

rule polishing_polypolish_insert_filter:
    """
    Exclude alignments based on insert size. 
    """
    input:
        SAM_1 = rules.polishing_read_mapping_1.output.SAM,
        SAM_2 = rules.polishing_read_mapping_2.output.SAM
    output:
        SAM = 'polish/short_reads/{assembly_type}/read_mapping/alignment_filtered_sam.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/read_mapping'
    threads: 8
    run:
        from camel.app.tools.polypolish.polypolishinsertfilter import PolypolishInsertFilter
        insert_filter = PolypolishInsertFilter()
        input_sam1 = snakemakeutils.load_object(Path(input.SAM_1))
        input_sam2 = snakemakeutils.load_object(Path(input.SAM_2))
        insert_filter.add_input_files({'SAM': [input_sam1[0], input_sam2[0]]})
        step = Step(rule_name=str(rule), tool=insert_filter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(insert_filter, output)

rule polishing_polypolish:
    """
    First polishing with polypolish.
    """
    input:
        SAM = rules.polishing_polypolish_insert_filter.output.SAM,
        FASTA = rules.polishing_samtools_index_polypolish.output.FASTA
    output:
        FASTA = 'polish/short_reads/{assembly_type}/polypolish/fasta.io',
        INFORMS = 'polish/short_reads/{assembly_type}/polypolish/informs.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/polypolish',
        polypolish_options = config.get('polishing', {}).get('polypolish', {})
    run:
        from camel.app.tools.polypolish.polypolish import Polypolish
        polypolish = Polypolish()
        snakemakeutils.add_pickle_inputs(polypolish, input)
        polypolish.update_parameters(**params.polypolish_options)
        step = Step(rule_name=str(rule), tool=polypolish, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(polypolish, output)

rule polishing_samtools_index_pypolca:
    """
    Creates a samtools index for the polypolish assembly.
    """
    input:
        FASTA = rules.polishing_polypolish.output.FASTA
    output:
        FASTA = 'polish/short_reads/{assembly_type}/pypolca/fasta-index.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/pypolca'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex()
        snakemakeutils.add_pickle_inputs(samtools, input)
        step = Step(rule_name=str(rule), tool=samtools, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools, output)

rule polishing_pypolca:
    """
    Then polishing with Pypolca.
    """
    input:
        FQ_dict = 'fq_dict.io',
        FASTA = rules.polishing_samtools_index_pypolca.output.FASTA
    output:
        FASTA = 'polish/short_reads/{assembly_type}/pypolca/fasta.io',
        INFORMS = 'polish/short_reads/{assembly_type}/pypolca/informs.io'
    params:
        dir_ = lambda wildcards: f'polish/short_reads/{wildcards.assembly_type}/pypolca',
        pypolca_options = config.get('polishing', {}).get('pypolca', {})
    threads: 8
    run:
        from camel.app.tools.pypolca.pypolca import Pypolca
        pypolca = Pypolca()
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict),'illumina')
        pypolca.add_input_files({'FASTQ_PE': fq_in.pe})
        snakemakeutils.add_pickle_input(pypolca, 'FASTA', Path(input.FASTA))
        pypolca.update_parameters(**params.pypolca_options, threads=threads)
        step = Step(rule_name=str(rule), tool=pypolca, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(pypolca, output)
