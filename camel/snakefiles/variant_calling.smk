from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import variant_calling, variant_filtering, read_simulation


rule variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        INDEX_GENOME_PREFIX = 'variant_calling/reference/genome_prefix.io',
        FASTA = 'variant_calling/reference/fasta.io',
        INFORMS = 'variant_calling/reference/informs.io'
    params:
        ref_fasta = config['reference'].get('fasta'),
        ref_url = config['reference'].get('url'),
        ref_name = config['reference'].get('name')
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.core.io.tooliofile import ToolIOFile
        snakemakeutils.dump_object([ToolIOValue(params.ref_fasta)], Path(output.INDEX_GENOME_PREFIX))
        snakemakeutils.dump_object([ToolIOFile(Path(params.ref_fasta))], Path(output.FASTA))
        snakemakeutils.dump_object({'name': params.ref_name, 'url': params.ref_url}, Path(output.INFORMS))

rule variant_calling_map_reads_illumina:
    """
    Maps the trimmed illumina reads to the reference sequence.
    """
    input:
        IO = variant_calling.get_mapping_fq_input(config),
        INDEX_GENOME_PREFIX = rules.variant_calling_prep_reference.output.INDEX_GENOME_PREFIX
    output:
        BAM = 'variant_calling/read_mapping/illumina/bam.io',
        INFORMS = 'variant_calling/read_mapping/illumina/informs.io'
    params:
        dir_ = 'variant_calling/read_mapping/illumina',
        input_type = config['input_type']
    threads: 4
    priority: 1
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map()
        if params.input_type != 'fasta':
            key_reads = 'PE' if params.input_type == 'illumina' else 'SE'
            bowtie2_map.add_input_files(snakepipelineutils.extract_fq_input(
                Path(input.IO), key_se='FASTQ_SE', drop_empty=True, read_type=key_reads))
        else:
            snakemakeutils.add_pickle_input(bowtie2_map, 'FASTQ_PE', Path(input.IO))
        bowtie2_map.update_parameters(threads=threads)
        snakemakeutils.add_pickle_input(bowtie2_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort], Path(params.dir_).absolute())

        # Save output
        snakemakeutils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        snakemakeutils.dump_object(bowtie2_map.informs, Path(output.INFORMS))

rule variant_calling_map_reads_ont:
    """
    Maps the trimmed nanopore reads to the reference sequence.
    """
    input:
        IO = variant_calling.get_mapping_fq_input(config),
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        BAM = 'variant_calling/read_mapping/ont/bam.io',
        INFORMS = 'variant_calling/read_mapping/ont/informs.io'
    params:
        dir_ = 'variant_calling/read_mapping/ont',
        input_type = config['input_type']
    threads: 4
    priority: 1
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Minimap2
        minimap2_map = Minimap2Mapping()
        snakemakeutils.add_pickle_input(minimap2_map, 'FASTA', Path(input.FASTA))
        minimap2_map.add_input_files(snakepipelineutils.extract_fq_input(Path(input.IO), key_se='FASTQ', read_type='SE'))
        minimap2_map.update_parameters(threads=threads)

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([minimap2_map, samtools_view, samtools_sort], Path(params.dir_).absolute())

        # Save output
        snakemakeutils.dump_tool_output(samtools_sort, 'BAM', Path(output.BAM))
        snakemakeutils.dump_object(minimap2_map.informs, Path(output.INFORMS))

rule variant_calling_generate_dummy_bam:
    """
    Generates a dummy bam.
    This rule is executed when the input type is FASTA with VCF.
    """
    output:
        BAM = 'variant_calling/dummy_bam/bam.io'
    params:
        dir_ = 'variant_calling/dummy_bam'
    run:
        import pysam

        header = {
            'HD': {'VN': '1.6'},
            'SQ': []  # No sequences
        }
        filename = Path(params.dir_, 'empty.bam')
        with pysam.AlignmentFile(str(filename), 'wb', header=header) as _:
            pass
        snakemakeutils.dump_object([ToolIOFile(filename)], Path(output.BAM))

rule variant_calling_calculate_mapping_rate:
    """
    Calculates the mapping rate of the input reads.
    """
    input:
        BAM = variant_calling.get_bam(config)
    output:
        INFORMS = variant_calling.OUTPUT_MAPPING_RATE_INFORMS
    params:
        dir_ = 'variant_calling/rate'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        samtools_flag = SamtoolsFlagstat()
        step = Step(rule_name=str(rule), tool=samtools_flag, dir_=Path(params.dir_))
        snakemakeutils.add_pickle_inputs(samtools_flag, input)
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_flag, output)

rule variant_calling_calculate_depth:
    """
    Calculates the median depth of the alignment.
    """
    input:
        BAM = variant_calling.get_bam(config)
    output:
        TSV = variant_calling.OUTPUT_DEPTH_TSV,
        INFORMS = variant_calling.OUTPUT_DEPTH_INFORMS
    params:
        dir_ = 'variant_calling/depth'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth()
        step = Step(rule_name=str(rule), tool=samtools_depth, dir_=Path(params.dir_))
        snakemakeutils.add_pickle_inputs(samtools_depth, input)
        samtools_depth.update_parameters(output_all_positions_absolutely=True)
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_depth, output)

rule variant_calling_mpileup:
    """
    This step creates a multi-way pileup using samtools.
    """
    input:
        FASTA = rules.variant_calling_prep_reference.output.FASTA,
        BAM = variant_calling.get_bam(config)
    output:
        VCF_GZ = 'variant_calling/mpileup/vcf_gz.io',
        INFORMS = 'variant_calling/mpileup/informs.io'
    priority: 1
    params:
        dir_ = 'variant_calling/mpileup',
        count_orphans = config['variant_calling'].get('count_orphans', True),
        min_mapping_quality = config['variant_calling'].get('minimal_mq'),
        min_base_quality = config['variant_calling'].get('minimal_bq'),
        disable_baq = config['variant_calling'].get('disable_baq'),
        input_type = config['input_type']
    threads: 1
    run:
        from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
        bcftools_mpileup = BcftoolsMpileup()
        snakemakeutils.add_pickle_inputs(bcftools_mpileup, input)
        step = Step(rule_name=str(rule), tool=bcftools_mpileup, dir_=Path(params.dir_))
        bcftools_mpileup.update_parameters(output_type='z')
        if params.input_type == 'ont':
            bcftools_mpileup.update_parameters(config='ont')
        if params.count_orphans is not None:
            bcftools_mpileup.update_parameters(count_orphans=params.count_orphans)
        if params.min_mapping_quality is not None:
            bcftools_mpileup.update_parameters(min_mapping_quality=params.min_mapping_quality)
        if params.min_base_quality is not None:
            bcftools_mpileup.update_parameters(min_base_quality=params.min_base_quality)
        if params.disable_baq is not None:
            bcftools_mpileup.update_parameters(disable_baq=params.disable_baq)
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_mpileup, output)

rule variant_calling_bcftools_call:
    """
    Variant calling using bcftools. 
    """
    input:
        VCF_GZ = rules.variant_calling_mpileup.output.VCF_GZ
    output:
        VCF_GZ = 'variant_calling/calling/vcf_gz.io',
        INFORMS = 'variant_calling/calling/informs.io'
    params:
        dir_ = 'variant_calling/calling',
        ploidy = config['variant_calling'].get('ploidy', 1),
        calling_method = config['variant_calling'].get('calling_method'),
        skip_variants = config['variant_calling'].get('skip_variants'),
        variants_only = config['variant_calling'].get('variants_only'),
        mutation_rate = config['variant_calling'].get('mutation_rate')
    run:
        from camel.app.tools.bcftools.bcftoolscall import BcftoolsCall
        variant_caller = BcftoolsCall()
        snakemakeutils.add_pickle_inputs(variant_caller, input)
        step = Step(rule_name=str(rule), tool=variant_caller, dir_=Path(str(params.dir_)))
        variant_caller.update_parameters(
            output_type='z',
            output_filename='variants.vcf.gz',
            variants_only=True,
            ploidy=params.ploidy
        )
        if params.calling_method is not None:
            variant_caller.update_parameters(calling_method=params.calling_method)
        if params.skip_variants is not None:
            variant_caller.update_parameters(skip_variants=params.skip_variants)
        if params.variants_only is not None:
            variant_caller.update_parameters(variants_only=params.variants_only)
        if params.mutation_rate is not None:
            variant_caller.update_parameters(mutation_rate=params.mutation_rate)
        step.run()
        snakemakeutils.dump_tool_outputs(variant_caller, output)

rule variant_calling_normalize_indels:
    """
    Normalizes indels.
    """
    input:
        VCF_GZ = rules.variant_calling_bcftools_call.output.VCF_GZ,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        VCF_GZ = 'variant_calling/norm/vcf_gz.io'
    params:
        dir_ = 'variant_calling/norm'
    run:
        from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm
        bcftools_norm = BcftoolsNorm()
        snakemakeutils.add_pickle_inputs(bcftools_norm, input)
        step = Step(rule_name=str(rule), tool=bcftools_norm, dir_=Path(str(params.dir_)))
        bcftools_norm.update_parameters(output_type='z')
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_norm, output)

rule variant_calling_index_vcf_gz:
    """
    Indexes the VCF file.
    """
    input:
        VCF_GZ = rules.variant_calling_normalize_indels.output.VCF_GZ
    output:
        VCF_GZ = variant_calling.OUTPUT_UNFILTERED_VCF_GZ
    params:
        dir_ = 'variant_calling/norm'
    run:
        from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
        indexer = BcftoolsIndex()
        snakemakeutils.add_pickle_inputs(indexer, input)
        step = Step(rule_name=str(rule), tool=indexer, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(indexer, output)

rule variant_calling_unzip_vcf:
    """
    Unzips the VCF file.
    """
    input:
        VCF_GZ = rules.variant_calling_index_vcf_gz.output.VCF_GZ
    output:
        VCF = variant_calling.OUTPUT_UNFILTERED_VCF
    params:
        dir_ = 'variant_calling/unzip_vcf',
        sample_name = config['sample_name']
    run:
        from camel.app.core.utils import fileutils
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView()
        snakemakeutils.add_pickle_inputs(bcftools_view, input)
        step = Step(rule_name=str(rule), tool=bcftools_view, dir_=Path(params.dir_))
        output_filename = f'variants-{fileutils.make_valid(params.sample_name)}.vcf'
        bcftools_view.update_parameters(output_type='v', output_filename=output_filename)
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_view, output)

rule variant_calling_zip_input_vcf:
    """
    Zips the input VCF file to be used by certain assays (e.g. amr detection) in case FASTA/VCF is used as input.
    """
    input:
        VCF = 'input/vcf.io'
    output:
        VCF_GZ = 'variant_calling/gzip/vcf_gz.io'
    params:
        running_dir = 'variant_calling/gzip'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView()
        snakemakeutils.add_pickle_inputs(bcftools_view, input)
        step = Step(rule_name=str(rule), tool=bcftools_view, dir_=Path(params.running_dir))
        bcftools_view.update_parameters(output_type='z')
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_view, output)

rule variant_calling_create_consensus:
    """
    Creates the consensus sequence by applying the detected variants to the reference genome.
    """
    input:
        VCF_GZ = rules.variant_calling_index_vcf_gz.output.VCF_GZ,
        FASTA = rules.variant_calling_prep_reference.output.FASTA
    output:
        FASTA = variant_calling.OUTPUT_CONSENSUS
    params:
        dir_ = 'variant_calling/consensus',
        output_filename = 'bcftools_consensus.fasta'
    run:
        from camel.app.tools.bcftools.bcftoolsconsensus import BcftoolsConsensus
        bcftools_consensus = BcftoolsConsensus()
        snakemakeutils.add_pickle_inputs(bcftools_consensus, input)
        bcftools_consensus.update_parameters(output_filename=params.output_filename)
        step = Step(rule_name=str(rule), tool=bcftools_consensus, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_consensus, output)

rule variant_calling_report:
    """
    Creates a report for the variant calling.
    """
    input:
        VCF = rules.variant_calling_unzip_vcf.output.VCF,
        VCF_filt = variant_filtering.OUTPUT_VCF,
        VCF_filt_regions = Path('variant_filtering', '06-regions', 'vcf.io') if variant_filtering.get_filtering_param(config, 'region', 'bed_file') is not None else [],
        BAM = variant_calling.get_bam(config),
        INFORMS_reference = rules.variant_calling_prep_reference.output.INFORMS,
        INFORMS_mapping = variant_calling.get_mapping_informs(config),
        INFORMS_calling = rules.variant_calling_bcftools_call.output.INFORMS,
        INFORMS_depth = rules.variant_calling_calculate_depth.output.INFORMS,
        INFORMS_map_rate = rules.variant_calling_calculate_mapping_rate.output.INFORMS,
        JSON = variant_filtering.OUTPUT_STATS
    output:
        VAL_HTML = variant_calling.OUTPUT_REPORT,
        INFORMS = 'variant_calling/report/informs.io'
    params:
        dir_ = 'variant_calling/report',
        regions_bed_file = variant_filtering.get_filtering_param(config, 'region', 'bed_file'),
        include_bam = config.get('variant_calling').get('report_include_bam', False),
        sample_name = config['sample_name'],
        input_type = config['input_type']
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.variant_calling.variantcallingreporter import VariantCallingReporter
        reporter = VariantCallingReporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        # noinspection PyUnresolvedReferences
        keys = [k for k in input.keys() if k != 'VCF_filt_regions']
        snakemakeutils.add_pickle_inputs(reporter, input, keys=keys)
        if params.regions_bed_file is not None:
            reporter.add_input_files({'BED': [ToolIOFile(Path(params.regions_bed_file))]})
            snakemakeutils.add_pickle_input(reporter, 'VCF_filt_regions', Path(input.VCF_filt_regions))
        reporter.update_parameters(export_bam='true' if params.include_bam else 'false')
        if params.input_type == 'fasta':
            reporter.update_parameters(pseudo_reads=True)
        reporter.add_input_files({'VAL_Sample': [ToolIOValue(params.sample_name)]})
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule variant_calling_dump_summary_info:
    """
    Dumps the summary information from the variant calling workflow.
    """
    input:
        INFORMS_map_rate = rules.variant_calling_calculate_mapping_rate.output.INFORMS,
        INFORMS_depth = rules.variant_calling_calculate_depth.output.INFORMS
    output:
        FILE = 'variant_calling/summary/summary_out.{ext}' # variant_calling.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs_map_rate = snakemakeutils.load_object(Path(input.INFORMS_map_rate))
        informs_depth = snakemakeutils.load_object(Path(input.INFORMS_depth))
        data_summary = [
            ('vc-mapping_rate', informs_map_rate['mapping_rate']),
            ('vc-median_depth', informs_depth['median_depth'])
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'variant_calling')

rule variant_calling_collect_command_informs:
    """
    This rule is used to collect the commands that were used.
    """
    input:
        INFORMS_read_simulation = read_simulation.OUTPUT_INFORMS if config['input_type'] == 'fasta' else [],
        INFORMS_mapping = variant_calling.get_mapping_informs(config),
        INFORMS_mpileup = rules.variant_calling_mpileup.output.INFORMS,
        INFORMS_calling = rules.variant_calling_bcftools_call.output.INFORMS
    output:
        INFORMS_ALL = 'variant_calling/informs_all.io', # variant_calling.OUTPUT_INFORMS_ALL,
        INFORMS = 'variant_calling/read_mapping/informs.io' # variant_calling.OUTPUT_MAPPING_INFORMS
    run:
        all_informs = []
        for io_file in input:
            informs = snakemakeutils.load_object(Path(io_file))
            informs['_tag'] = 'Variant calling'
            all_informs.append(informs)
        snakemakeutils.dump_object(all_informs, Path(output.INFORMS_ALL))
        # dump the mapping informs to file (needed in assembly.py)
        informs_map = snakemakeutils.load_object(Path(input.INFORMS_mapping))
        snakemakeutils.dump_object(informs_map, Path(output.INFORMS))

rule variant_calling_report_empty:
    """
    Creates an empty variant calling report when a VCF file is given as input.
    """
    output:
        HTML = 'variant_calling/report/html-empty.iob' # variant_calling.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Variant calling', Path(output.HTML))
