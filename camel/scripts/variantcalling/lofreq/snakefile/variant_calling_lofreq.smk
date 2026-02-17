from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.tools.lofreq.lofreqreporter import LofreqReporter
from camel.snakefiles import trimming_illumina, trimming
from camel.scripts.variantcalling.lofreq.snakefile import variant_calling_lofreq
from camel.app.core.reports import reportutils

include: trimming_illumina.SNAKEFILE

rule all:
    input:
        VCF='variant_calling/vcf.io',
        HTML='variant_calling/report/html.iob'

rule link_fastq_to_trimming_input:
    """
    Creates the FASTQ input for the human read scrubbing step.
    """
    output:
        FASTQ_PE=trimming_illumina.INPUT_FASTQ,
    params:
        input_dict=config['input']
    run:
        from camel.app.core.io.tooliofile import ToolIOFile

        snakemakeutils.dump_object([ToolIOFile(Path(x['path'])) for x in
                                    params.input_dict['fastq_pe']],Path(output.FASTQ_PE))

rule variant_calling_prep_reference:
    """
    Converts the reference to a Snakemake / CAMEL compatible format. Creates the reference metadata.
    """
    output:
        INDEX_GENOME_PREFIX='variant_calling/reference/genome_prefix.io',
        FASTA='variant_calling/reference/fasta.io',
        INFORMS='variant_calling/reference/informs.io'
    params:
        ref_fasta=config['reference'].get('fasta'),
        ref_url=config['reference'].get('url'),
        ref_name=config['reference'].get('name')
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.core.io.tooliofile import ToolIOFile

        snakemakeutils.dump_object([ToolIOValue(params.ref_fasta)],Path(output.INDEX_GENOME_PREFIX))
        snakemakeutils.dump_object([ToolIOFile(Path(params.ref_fasta))],Path(output.FASTA))
        snakemakeutils.dump_object({'name': params.ref_name, 'url': params.ref_url},Path(output.INFORMS))

rule variant_calling_map_reads_illumina:
    """
    Maps the trimmed illumina reads to the reference sequence.
    """
    input:
        IO=trimming_illumina.select_fastq_output(config),
        INDEX_GENOME_PREFIX=rules.variant_calling_prep_reference.output.INDEX_GENOME_PREFIX
    output:
        BAM='variant_calling/read_mapping/illumina/bam.io',
        INFORMS='variant_calling/read_mapping/illumina/informs.io'
    params:
        dir_='variant_calling/read_mapping/illumina',
        input_type=config['input']['type']
    threads: 4
    priority: 1
    run:
        from camel.app.core.piping import pipeutils
        from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        # Bowtie 2
        bowtie2_map = Bowtie2Map()
        snakemakeutils.add_pickle_input(bowtie2_map,'FASTQ_PE',Path(input.IO))
        bowtie2_map.update_parameters(threads=threads)
        snakemakeutils.add_pickle_input(bowtie2_map,'INDEX_GENOME_PREFIX',Path(input.INDEX_GENOME_PREFIX))

        # Initialize tools
        samtools_view = SamtoolsView()
        samtools_sort = SamtoolsSort()
        samtools_sort.update_parameters(threads=threads)
        pipeutils.run_as_pipe([bowtie2_map, samtools_view, samtools_sort],Path(params.dir_).absolute())

        # Save output
        snakemakeutils.dump_tool_output(samtools_sort,'BAM',Path(output.BAM))
        snakemakeutils.dump_object(bowtie2_map.informs,Path(output.INFORMS))

rule lofreq_indel_qualities:
    input:
        BAM=rules.variant_calling_map_reads_illumina.output.BAM,
        FASTA=rules.variant_calling_prep_reference.output.FASTA
    output:
        BAM='variant_calling/read_mapping/illumina/bam-indelqual.io',
        INFORMS='variant_calling/read_mapping/illumina/indelqual_informs.io'
    params:
        dir_='variant_calling/read_mapping/illumina'
    run:
        from camel.app.tools.lofreq.lofreqindelqual import LofreqIndelqual

        lofreq_indelqual = LofreqIndelqual()
        snakemakeutils.add_pickle_inputs(lofreq_indelqual,input)
        step = Step(rule_name=str(rule),tool=lofreq_indelqual,dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(lofreq_indelqual,output)

rule variant_calling_with_lofreq:
    input:
        BAM=rules.lofreq_indel_qualities.output.BAM if config['variant_calling']['call_indels'] is True
        else rules.variant_calling_map_reads_illumina.output.BAM,
        FASTA=rules.variant_calling_prep_reference.output.FASTA
    output:
        VCF='variant_calling/vcf.io',
        INFORMS='variant_calling/informs.io',
    params:
        dir_=lambda wildcards: 'variant_calling/'
    threads: 8
    run:
        from camel.app.tools.lofreq.lofreqcall import LofreqCall

        lofreq_call = LofreqCall()
        snakemakeutils.add_pickle_inputs(lofreq_call,input)
        step = Step(rule_name=str(rule),tool=lofreq_call,dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(lofreq_call,output)

rule lofreq_reporter:
    input:
        VCF=rules.variant_calling_with_lofreq.output.VCF,
        INFORMS_reference=rules.variant_calling_prep_reference.output.INFORMS
    output:
        VAL_HTML=variant_calling_lofreq.OUTPUT_REPORT_LOFREQ,
        INFORMS='variant_calling/report/informs.io'
    params:
        dir_='variant_calling/report',
        sample_name=config['input']['sample_name']
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue

        reporter = LofreqReporter()
        step = Step(rule_name=str(rule),tool=reporter,dir_=Path(params.dir_))
        snakemakeutils.add_pickle_inputs(reporter,input)
        reporter.add_input_files({'VAL_Sample': [ToolIOValue(params.sample_name)]})
        step.run()
        snakemakeutils.dump_tool_outputs(reporter,output)

rule report_create_commands_section:
    """
    Creates the section with the commands.
    """
    input:
        INFORMS_trimming=trimming.get_command_informs(config),
        INFORMS_indelqual=variant_calling_lofreq.OUTPUT_INDELQUAL_INFORMS if config['variant_calling']['call_indels'] is True else [],
        INFORMS_lofreq=variant_calling_lofreq.OUTPUT_LOFREQ_INFORMS
    output:
        HTML='report/html-commands.iob'
    params:
        dir_=config['working_dir']
    run:
        from camel.app.scriptutils.basepipe import basepipeutils

        basepipeutils.export_command_section(input,Path(output.HTML),params.dir_)

rule generate_report:
    input:
        report_trimming=trimming.get_reports(config),
        report_lofreq=variant_calling_lofreq.OUTPUT_REPORT_LOFREQ,
        report_commands=rules.report_create_commands_section.output.HTML
    output:
        HTML=config['output']['html']
    params:
        output_dir=config['output']['dir'],
        html=config['output']['html'],
        name=config['script_info']['name'],
        version=config['script_info']['version'],
        input_dict=config['input']
    run:
        import datetime
        from camel.app.scriptutils.basepipe import basepipeutils
        from camel.app.scriptutils.basescript.scriptinput import ScriptInput

        # Initialize report
        report = reportutils.init_report(
            path_out=Path(params.html),
            key=params.name,
            title=params.name,
        )
        report.add_html_object(reportutils.create_overview_section(
            version=params.version,
            dataset_name=params.input_dict['name'],
            input_file_str=params.input_dict['input_str'],
            date=datetime.datetime.now()
        ))

        report_structure = []

        report_structure.append(('Trimming', 'trimming', [Path(input.report_trimming[0])]))
        report_structure.append(('Variant Calling', 'var_calling', [Path(input.report_lofreq)]))
        report_structure.append(('Commands', 'commands', [Path(input.report_commands)]))
        snakepipelineutils.add_report_content(report,report_structure)
        report.save()
