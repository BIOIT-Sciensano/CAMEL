from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import polish_assembly_long


rule medaka_polishing_map_ont_reads:
    """
    Maps the ONT reads to the Flye assembly.
    """
    input:
        FQ = 'fq_dict.io',
        FASTA = lambda wildcards: polish_assembly_long.INPUT_ASSEMBLY_FASTA.format(assembly_type=wildcards.assembly_type)
    output:
        BAM = 'polish/long_reads/{assembly_type}/minimap2/bam.io'
    params:
        dir_ = lambda wildcards: f'polish/long_reads/{wildcards.assembly_type}/minimap2'
    threads: 8
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Minimap2
        minimap2 = Minimap2Mapping()
        snakemakeutils.add_pickle_input(minimap2, 'FASTA', Path(str(input.FASTA)))
        minimap2.add_input_files(snakepipelineutils.extract_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([minimap2, samtools_view, samtools_sort], Path(str(params.dir_)))

        # Export output
        snakemakeutils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))

rule medaka_polishing_index_bam_file:
    input:
        BAM = rules.medaka_polishing_map_ont_reads.output.BAM
    output:
        BAM = 'polish/long_reads/{assembly_type}/minimap2/bam-index.io'
    params:
        dir_ = lambda wildcards: f'polish/long_reads/{wildcards.assembly_type}/minimap2'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex()
        snakemakeutils.add_pickle_input(samtools_index, 'BAM', Path(input.BAM))
        samtools_index.update_parameters(generate_bai_index=True)
        step = Step(rule_name=str(rule), tool=samtools_index, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_index, output)

rule medaka_polishing_medaka_inference:
    """
    Runs the medaka inference algorithm.
    """
    input:
        BAM = rules.medaka_polishing_index_bam_file.output.BAM
    output:
        HDF = 'polish/long_reads/{assembly_type}/inference/raw_hdf.io',
        INFORMS = 'polish/long_reads/{assembly_type}/inference/commands-inference.io'
    params:
        dir_ =  lambda wildcards: f'polish/long_reads/{wildcards.assembly_type}/inference',
        medaka_options = config.get('polishing', {}).get('medaka', {}).get('inference', {})
    threads: 8
    run:
        from camel.app.tools.medaka.medakainference import MedakaInference
        medaka = MedakaInference()
        snakemakeutils.add_pickle_input(medaka, 'BAM', Path(input.BAM))
        medaka.update_parameters(**params.medaka_options)
        medaka.update_parameters(threads=threads)
        step = Step(rule_name=str(rule), tool=medaka, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(medaka, output)

rule medaka_polishing_medaka_sequence:
    """
    Runs the medaka sequence algorithm.
    """
    input:
        HDF = rules.medaka_polishing_medaka_inference.output.HDF,
        FASTA = 'polish/long_reads/{assembly_type}/input/fasta.io'
    output:
        FASTA = 'polish/long_reads/{assembly_type}/sequence/fasta.io',
        INFORMS = 'polish/long_reads/{assembly_type}/sequence/commands-sequence.io'
    params:
        dir_ = lambda wildcards: f'polish/long_reads/{wildcards.assembly_type}/sequence',
        medaka_options = config.get('polishing', {}).get('medaka', {}).get('sequence', {})
    threads: 8
    run:
        from camel.app.tools.medaka.medakasequence import MedakaSequence
        medaka = MedakaSequence()
        snakemakeutils.add_pickle_inputs(medaka, input)
        medaka.update_parameters(**params.medaka_options, threads=threads)
        step = Step(rule_name=str(rule), tool=medaka, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(medaka, output)

rule medaka_polishing_empty_report:
    """
    Creates an empty report for the gene detection when plasmidSPAdes assembly fails.
    """
    output:
        VAL_HTML = 'polish/long_reads/{assembly_type}/report/html-empty.iob'
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Medaka polishing', Path(output.VAL_HTML))
