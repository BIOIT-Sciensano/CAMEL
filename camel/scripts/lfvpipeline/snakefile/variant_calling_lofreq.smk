from pathlib import Path

from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.core.snakemake.step import Step
from camel.app.tools.lofreq.lofreqreporter import LofreqReporter
from camel.scripts.lfvpipeline.snakefile import variant_calling_lofreq
from camel.snakefiles import trimming_illumina, trimming, variant_calling

include: trimming_illumina.SNAKEFILE
include: variant_calling.SNAKEFILE

rule link_fastq_to_trimming_input:
    """
    Creates the FASTQ input for the human read scrubbing step.
    """
    output:
        FASTQ_PE = trimming_illumina.INPUT_FASTQ
    params:
        input_dict = config['input']
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        snakemakeutils.dump_object([
            ToolIOFile(Path(x['path'])) for x in params.input_dict['fastq_pe']], Path(output.FASTQ_PE))

rule fai_index:
    """
    Creates a bowtie2 index for the reference genome.
    """
    input:
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        FASTA = 'variant_calling/reference/genome_prefix_index_fai.io'
    params:
        dir_ = 'variant_calling/reference'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex

        samtools_index = SamtoolsFastaIndex()
        step = Step(rule_name=str(rule), tool=samtools_index, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(samtools_index, input)
        step.run()
        snakemakeutils.dump_io_outputs(samtools_index, output)

rule create_sequence_dictionary:
    """
    Creates a sequence dictionary for the reference genome.
    """
    input:
        FASTA_REF = rules.fai_index.output.FASTA
    output:
        FASTA_REF = 'variant_calling/reference/fasta_sequence_dictionary.io'
    params:
        dir_ = 'variant_calling/reference'
    run:
        from camel.app.tools.picard.createsequencedictionary import CreateSequenceDictionary

        create_dictionary = CreateSequenceDictionary()
        step = Step(rule_name=str(rule), tool=create_dictionary, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(create_dictionary, input)
        step.run()
        snakemakeutils.dump_io_outputs(create_dictionary, output)

rule bt2_index:
    """
    Creates a bowtie2 index for the reference genome.
    """
    input:
        FASTA_REF = rules.fai_index.output.FASTA
    output:
        INDEX_GENOME_PREFIX = 'variant_calling/reference/genome_prefix_index.io'
    params:
        dir_ = 'variant_calling/reference'
    run:
        from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index

        bowtie2_index = Bowtie2Index()
        step = Step(rule_name=str(rule), tool=bowtie2_index, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(bowtie2_index, input)
        step.run()
        snakemakeutils.dump_io_outputs(bowtie2_index, output)

rule variant_calling_map_reads_illumina_lofreq:
    """
    Maps the trimmed illumina reads to the reference sequence.
    """
    input:
        IO = trimming_illumina.select_fastq_output(config),
        INDEX_GENOME_PREFIX = rules.bt2_index.output.INDEX_GENOME_PREFIX
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam.io',
        INFORMS = 'variant_calling/read_mapping/illumina/informs.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina',
        input_type = config['input']['type']
    threads: 4
    priority: 1
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map()
        snakemakeutils.add_io_input(bowtie2_map, 'FASTQ_PE', Path(input.IO))
        bowtie2_map.update_parameters(threads=threads)
        snakemakeutils.add_io_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.dir_).absolute())

        # Save output
        snakemakeutils.dump_io_output(samtools_sort, 'BAM', Path(output.BAM))
        snakemakeutils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule picard_add_readgroups:
    """
    Add readgroup information to the BAM files. Necessary for the IndelRealigner
    """
    input:
        BAM = rules.variant_calling_map_reads_illumina_lofreq.output.BAM
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam-rg.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina'
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups

        add_rg = AddOrReplaceReadGroups()
        snakemakeutils.add_io_inputs(add_rg, input)
        step = Step(rule_name=str(rule), tool=add_rg, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(add_rg, output)

rule bam_index:
    """
    Add readgroup information to the BAM files.
    """
    input:
        BAM = rules.picard_add_readgroups.output.BAM
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam-rg-index.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex

        samtools_index = SamtoolsIndex()
        snakemakeutils.add_io_inputs(samtools_index, input)
        step = Step(rule_name=str(rule), tool=samtools_index, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(samtools_index, output)

rule gatk_realigner_target_creator:
    input:
        BAM = rules.bam_index.output.BAM,
        FASTA_REF = rules.create_sequence_dictionary.output.FASTA_REF
    output:
        TXT_realign_intervals = 'variant_calling/read_mapping/illumina/txt-realign-intervals.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina'
    run:
        from camel.app.tools.gatk.gatkrealignertargetcreator import GATKRealignerTargetCreator
        gatk_realigner = GATKRealignerTargetCreator()
        snakemakeutils.add_io_inputs(gatk_realigner, input)
        step = Step(rule_name=str(rule), tool=gatk_realigner, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(gatk_realigner, output)

rule gatk_indel_realigner:
    input:
        BAM = rules.bam_index.output.BAM,
        FASTA_REF = rules.create_sequence_dictionary.output.FASTA_REF,
        TXT_realign_intervals = rules.gatk_realigner_target_creator.output.TXT_realign_intervals
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam_realigned/bam.io',
        INFORMS = 'variant_calling/read_mapping/illumina/bam_realigned/informs.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina/bam_realigned'
    run:
        from camel.app.tools.gatk.gatkindelrealigner import GATKIndelRealigner
        gatk_realigner = GATKIndelRealigner()
        snakemakeutils.add_io_inputs(gatk_realigner,input)
        step = Step(rule_name=str(rule),tool=gatk_realigner,dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(gatk_realigner,output)

rule lofreq_indel_qualities:
    """
    Inserts indel qualities into the BAM file, using LoFreq indelqual, necessary for calling indels.
    """
    input:
        BAM = rules.gatk_indel_realigner.output.BAM,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam-indelqual.io',
        INFORMS = 'variant_calling/read_mapping/illumina/indelqual_informs.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina'
    run:
        from camel.app.tools.lofreq.lofreqindelqual import LofreqIndelqual

        lofreq_indelqual = LofreqIndelqual()
        snakemakeutils.add_io_inputs(lofreq_indelqual, input)
        step = Step(rule_name=str(rule), tool=lofreq_indelqual, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(lofreq_indelqual, output)

rule variant_calling_with_lofreq:
    """
    Executes the variant calling with LoFreq call.
    """
    input:
        BAM = rules.lofreq_indel_qualities.output.BAM if config['variant_calling']['call_indels'] is True
        else rules.gatk_indel_realigner.output.BAM,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        VCF = 'variant_calling/vcf.io',
        INFORMS = 'variant_calling/informs.io'
    params:
        dir_ = 'variant_calling',
        call_indels = config.get('variant_calling', {}).get('call_indels', None)
    threads: 8
    run:
        from camel.app.tools.lofreq.lofreqcall import LofreqCall

        lofreq_call = LofreqCall()
        snakemakeutils.add_io_inputs(lofreq_call, input)
        lofreq_call.update_parameters(call_indels=True if params.call_indels else False)
        step = Step(rule_name=str(rule), tool=lofreq_call, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(lofreq_call, output)

rule annotate_variants_csq:
    """
    Determines the consequence of the mutation using bcftools csq.
    """
    input:
        VCF = rules.variant_calling_with_lofreq.output.VCF,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        VCF = 'variant_calling/csq/vcf_csq.io',
        INFORMS = 'variant_calling/csq/csq_informs.io'
    params:
        dir_ = 'variant_calling/csq/',
        gff = config.get('reference', {}).get('gff', None)
    run:
        from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
        from camelcore.app.io.tooliofile import ToolIOFile

        csq = BcftoolsCsq()
        snakemakeutils.add_io_inputs(csq, input)
        if params.gff:
            csq.add_input_files({'GFF': [ToolIOFile(Path(params.gff))]})
        step = Step(rule_name=str(rule), tool=csq, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(csq, output)

rule extract_variants_effect_from_vcf:
    """
    Extracts all variants and the consequence (if applied) and stores them in a TSV.
    """
    input:
        VCF = rules.annotate_variants_csq.output.VCF if config['variant_calling'].get('csq', False) is True
        else rules.variant_calling_with_lofreq.output.VCF
    output:
        TSV = 'variant_calling/variants_list/tsv.iob'
    params:
        dir_ = 'variant_calling/variants_list'
    run:
        from camel.app.tools.pipelines.variant_calling.extractvariantsandeffectfromvcf import ExtractVariantsAndEffectFromVCF
        extract_variants = ExtractVariantsAndEffectFromVCF()
        step = Step(rule_name=str(rule), tool=extract_variants, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(extract_variants, input)
        step.run()
        snakemakeutils.dump_io_outputs(extract_variants, output)

rule lofreq_reporter:
    input:
        VCF = rules.annotate_variants_csq.output.VCF if config['variant_calling'].get('csq', False) is True
        else rules.variant_calling_with_lofreq.output.VCF,
        INFORMS_reference = rules.variant_calling_prep_reference.output.INFORMS,
        INFORMS_mapping = variant_calling.get_mapping_informs(config),
        INFORMS_map_rate = rules.variant_calling_calculate_mapping_rate.output.INFORMS,
        INFORMS_depth = rules.variant_calling_calculate_depth.output.INFORMS,
        INFORMS_lofreq = rules.variant_calling_with_lofreq.output.INFORMS,
        BAM = 'variant_calling/read_mapping/illumina/bam.io',
        TSV_depth = rules.variant_calling_calculate_depth.output.TSV,
        TSV_list = rules.extract_variants_effect_from_vcf.output.TSV
    output:
        VAL_HTML = 'variant_calling/report/html-lofreq.iob',  # variant_calling_lofreq.OUTPUT_REPORT_LOFREQ
        INFORMS = 'variant_calling/report/informs.io'
    params:
        dir_ = 'variant_calling/report',
        sample_name = config['input']['sample_name'],
        include_bam = config.get('variant_calling', {}).get('report_include_bam', False),
        min_af = config.get('variant_calling', {}).get('min_af', 0),
        csq = config.get('variant_calling', {}).get('csq', False)
    run:
        reporter = LofreqReporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        snakemakeutils.add_io_inputs(reporter, input)
        reporter.update_parameters(
            export_bam='true' if params.include_bam else 'false',
            min_af=params.min_af,
            sample_name=params.sample_name)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule report_create_commands_section:
    """
    Creates the section with the commands.
    """
    input:
        INFORMS_trimming = trimming.get_command_informs(config),
        INFORMS_indelqual = variant_calling_lofreq.OUTPUT_INDELQUAL_INFORMS if
        config.get('variant_calling', {}.get('call_indels', False)) is True else [],
        INFORMS_lofreq = variant_calling_lofreq.OUTPUT_LOFREQ_INFORMS
    output:
        HTML = 'report/html-commands.iob'
    params:
        dir_ = config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils

        basepipeutils.export_command_section(input, Path(output.HTML), params.dir_)

rule generate_report:
    input:
        report_trimming = trimming.get_reports(config),
        report_lofreq = variant_calling_lofreq.OUTPUT_REPORT_LOFREQ,
        report_commands = rules.report_create_commands_section.output.HTML
    output:
        HTML = config['output']['html']
    params:
        output_dir = config['output']['dir'],
        name = config['script_info']['name'],
        version = config['script_info']['version'],
        input_dict = config['input'],
        pipeline_info = config['script_info']
    run:
        from camelcore.app.utils import reportutils

        # Initialize report
        report = snakepipelineutils.init_pipeline_report(
            Path(output.HTML),Path(params.output_dir), params.pipeline_info)
        report.add_html_object(reportutils.create_overview_section(
            version=params.version,
            dataset_name=params.input_dict['name'],
            input_file_str=params.input_dict['input_str']
        ))
        report_structure = [
            ('Trimming', 'trimming', [Path(input.report_trimming[0])]),
            ('Variant Calling', 'var_calling', [Path(input.report_lofreq)]),
            ('Commands', 'commands', [Path(input.report_commands)])
        ]
        snakepipelineutils.add_report_content(report, report_structure)
        report.save()
