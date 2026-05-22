from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly, assembly_spades, assembly_flye, polish_assembly_short, polish_assembly_long


include: assembly_spades.SNAKEFILE
include: assembly_flye.SNAKEFILE
include: polish_assembly_long.SNAKEFILE
include: polish_assembly_short.SNAKEFILE


rule assembly_filter_contig_length:
    """
    Filters out the small contigs.
    """
    input:
        FASTA = assembly.get_fasta_raw(config)
    output:
        FASTA = 'assembly/filtering/fasta.io', # assembly.OUTPUT_FASTA
        INFORMS = 'assembly/filtering/informs.io' # assembly.OUTPUT_INFORMS_FILTERING
    params:
        dir_ = 'assembly/filtering',
        min_contig_length = config.get('assembly', {}).get('min_contig_length', 0)
    run:
        from camel.app.tools.seqtk.seqtkseq import SeqtkSeq
        seqtk = SeqtkSeq()
        snakemakeutils.add_io_inputs(seqtk, input)
        step = Step(rule_name=str(rule), tool=seqtk, dir_=Path(params.dir_))
        seqtk.update_parameters(output_filename='assembly_filtered.fasta', min_length=params.min_contig_length)
        step.run()
        snakemakeutils.dump_io_outputs(seqtk, output)

rule assembly_bt2_index:
    """
    Creates a bowtie2 index for the assembly.
    """
    input:
        FASTA_REF = rules.assembly_filter_contig_length.output.FASTA
    output:
        INDEX_GENOME_PREFIX = 'assembly/bowtie2/genome_prefix.io'
    params:
        dir_ = 'assembly/bowtie2'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
        bowtie2_index = Bowtie2Index()
        step = Step(rule_name=str(rule), tool=bowtie2_index, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(bowtie2_index, input)
        step.run()
        snakemakeutils.dump_io_outputs(bowtie2_index, output)

rule assembly_bt2_map:
    """
    Maps the reads against the assembled contigs.
    """
    input:
        IO = 'fq_dict.io',
        INDEX_GENOME_PREFIX = rules.assembly_bt2_index.output.INDEX_GENOME_PREFIX
    output:
        BAM = 'assembly/bowtie2/bam.io',
        INFORMS = 'assembly/bowtie2/informs.io' # assembly.OUTPUT_INFORMS_MAPPING
    params:
        dir_ = 'assembly/bowtie2',
        read_type = 'SE' if config.get('read_type') == 'iontorrent' else 'PE'
    threads: 8
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map()
        bowtie2_map.add_input_files(snakepipelineutils.extract_fq_input(
            Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=params.read_type))
        bowtie2_map.update_parameters(threads=threads)
        snakemakeutils.add_io_input(bowtie2_map,'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.dir_))

        # Save output
        snakemakeutils.dump_io_output(samtools_sort,'BAM', Path(output.BAM))
        bowtie2_map.informs['_tag'] = 'Coverage calculation'
        snakemakeutils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule assembly_minimap2_map:
    """
    Maps the reads to the Flye assembly using minimap2.
    """
    input:
        FQ = 'fq_dict.io',
        FASTA = rules.assembly_filter_contig_length.output.FASTA
    output:
        BAM = 'assembly/minimap2/bam.io',
        INFORMS = 'assembly/minimap2/informs.io' # assembly.OUTPUT_INFORMS_MAPPING
    params:
        dir_ = 'assembly/minimap2'
    threads: 8
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Minimap2
        minimap2 = Minimap2Mapping()
        snakemakeutils.add_io_input(minimap2,'FASTA', Path(input.FASTA))
        minimap2.add_input_files(snakepipelineutils.extract_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(rule_name=str(rule), tool=minimap2, dir_=Path(params.dir_))
        step.run()

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([minimap2, samtools_view, samtools_sort], Path(params.dir_))

        # Export output
        snakemakeutils.dump_io_output(samtools_sort,'BAM', Path(output.BAM))
        minimap2.informs['_tag'] = 'Coverage calculation'
        snakemakeutils.dump_object(minimap2.informs, Path(output.INFORMS))

rule assembly_samtools_depth:
    """
    Runs samtools depth on the BAM file of the reads mapped to the assembly.
    """
    input:
        BAM = 'assembly/{mapper}/bam.io'
    output:
        TSV = 'assembly/samtools_depth/{mapper}/tsv.io',
        INFORMS = 'assembly/samtools_depth/{mapper}/informs.io' # assembly.OUTPUT_INFORMS_DEPTH
    params:
        dir_ = lambda wildcards: f'assembly/samtools_depth/{wildcards.mapper}'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth()
        step = Step(rule_name=str(rule), tool=samtools_depth, dir_=Path(str(params.dir_)))
        snakemakeutils.add_io_inputs(samtools_depth, input)
        step.run()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        snakemakeutils.dump_io_outputs(samtools_depth, output)

rule assembly_samtools_flagstat:
    """
    Runs samtools flagstat to determine the mapping rate.
    """
    input:
        BAM = 'assembly/{mapper}/bam.io'
    output:
        INFORMS = 'assembly/samtools_flagstat/{mapper}/informs.io' # assembly.OUTPUT_INFORMS_MAPPING_RATE
    params:
        dir_ = 'assembly/samtools_flagstat'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        flagstat = SamtoolsFlagstat()
        snakemakeutils.add_io_inputs(flagstat, input)
        step = Step(rule_name=str(rule), tool=flagstat, dir_=Path(params.dir_))
        step.run()
        # Calculate mapping rate
        flagstat.informs['mapping_perc'] = 100 * flagstat.informs['mapped'][0] / flagstat.informs['total'][0]
        snakemakeutils.dump_io_outputs(flagstat, output)

#################################################################
# The rules below generate a basic QUAST report / summary file. #
# For a more in-depth QUAST report use the quast.smk snakefile  #
#################################################################
rule assembly_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.assembly_filter_contig_length.output.FASTA
    output:
        TSV = 'assembly/quast/tsv.io',
        INFORMS = 'assembly/quast/informs.io'
    params:
        dir_ = 'assembly/quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast()
        snakemakeutils.add_io_inputs(quast, input)
        step = Step(rule_name=str(rule), tool=quast, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(quast, output)

rule assembly_quast_extract_informs:
    """
    Extracts the information from the QUAST output file.
    """
    input:
        TSV = rules.assembly_quast.output.TSV,
        INFORMS_quast = rules.assembly_quast.output.INFORMS
    output:
        INFORMS = 'assembly/quast_extract/informs.io'
    params:
        dir_ = 'assembly/quast_extract'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor()
        snakemakeutils.add_io_inputs(quast_inform_extractor, input)
        step = Step(rule_name=str(rule), tool=quast_inform_extractor, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(quast_inform_extractor, output)

rule assembly_quast_report:
    """
    Creates a basic HTML report for the assembly.
    """
    input:
        FASTA_raw = assembly.get_fasta_raw(config),
        FASTA_filt = rules.assembly_filter_contig_length.output.FASTA,
        INFORMS_assembler = assembly.get_command_informs(config),
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        VAL_HTML = 'assembly/report/html.iob'
    params:
        dir_ = 'assembly/report',
        sample_name = config['input']['sample_name']
    run:
        from camel.app.tools.pipelines.assembly.htmlreporterassembly import HtmlReporterAssembly
        from camelcore.app.io.tooliovalue import ToolIOValue
        reporter = HtmlReporterAssembly()
        reporter.add_input_files({'SAMPLE_NAME': [ToolIOValue(params.sample_name)]})
        assembler_name = ', '.join(snakemakeutils.load_object(Path(x))['_name'] for x in input.INFORMS_assembler) \
            if input.INFORMS_assembler else 'n/a'
        reporter.add_input_informs({'assembler': {'_name': assembler_name}})
        snakemakeutils.add_io_inputs(reporter, input, excluded_keys=['INFORMS_assembler'])
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule assembly_dump_summary_info:
    """
    Dumps the summary information from the assembly pipeline.
    """
    input:
        INFORMS_quast = rules.assembly_quast_extract_informs.output.INFORMS
    output:
        TSV = assembly.OUTPUT_SUMMARY
    run:
        quast_informs = snakemakeutils.load_object(Path(input.INFORMS_quast))
        summary_data = [
            ('assembly_n50', quast_informs['contig']['N50']),
            ('assembly_nb_contigs', quast_informs['contig']['# contigs']),
            ('assembly_total_length', quast_informs['genome']['Total length'])
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in summary_data:
                handle.write(f'{key}\t{value}')
                handle.write('\n')
