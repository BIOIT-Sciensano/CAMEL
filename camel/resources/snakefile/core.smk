import shutil
from pathlib import Path

from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import trimming_illumina, trimming_ont, assembly, gene_detection, downsampling, core, \
    sequence_typing, amrfinder, mobsuite, assembly_flye, polish_assembly_short, polish_assembly_long, \
    human_read_scrubbing


rule core_link_fastq_scrubbing_input:
    """
    Creates the FASTQ input for the human read scrubbing step.
    """
    output:
        FASTQ_PE = human_read_scrubbing.INPUT_FASTQ.format(input_format='fastq_pe'),
        FASTQ_SE = human_read_scrubbing.INPUT_FASTQ.format(input_format='fastq_se')
    params:
        input_dict = config['input']
    run:
        from camel.app.io.tooliofile import ToolIOFile

        # FASTQ PE
        if 'fastq_pe' in params.input_dict:
            snakemakeutils.dump_object([ToolIOFile(Path(x['path'])) for x in params.input_dict['fastq_pe']], Path(output.FASTQ_PE))
        else:
            snakemakeutils.dump_object([], Path(output.FASTQ_PE))

        # FASTQ SE
        if 'fastq_se' in params.input_dict:
            snakemakeutils.dump_object([ToolIOFile(Path(params.input_dict['fastq_se'][0]['path']))], Path(output.FASTQ_SE))
        else:
            snakemakeutils.dump_object([], Path(output.FASTQ_SE))

rule core_link_fastq_pe_downsampling_input:
    """
    Links the human read scrubbing input or output to the input of the downsampling workflow.
    """
    input:
        FASTQ = core.get_fastq_input_downsampling(config, 'fastq_pe')
    output:
        FASTQ = downsampling.INPUT_FASTQ.format(read_key='fastq_pe')
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_link_fastq_se_downsampling_input:
    """
    Links the human read scrubbing input to the downsampling workflow.
    If human read scrubbing is disabled, the pipeline input is linked instead. 
    """
    input:
        FASTQ = core.get_fastq_input_downsampling(config, 'fastq_se')
    output:
        FASTQ = downsampling.INPUT_FASTQ.format(read_key='fastq_se')
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_link_trimmomatic_input:
    """
    Links the downsampling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = downsampling.OUTPUT_FASTQ.format(read_key='fastq_pe')
    output:
        FASTQ = trimming_illumina.INPUT_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_link_ont_input:
    """
    Links the downsampling output to the input of the trimmomatic workflow.  
    """
    input:
        FASTQ = downsampling.OUTPUT_FASTQ.format(read_key='fastq_se')
    output:
        FASTQ = trimming_ont.INPUT_ONT_FASTQ
    shell:
        """
        cp {input.FASTQ} {output.FASTQ};
        """

rule core_collect_trimmed_fastq_data:
    """
    This rule core_creates an IO object with the trimmed FASTQ files.
    Other workflows such as 'Kraken' or 'Assembly' rely on this dictionary to get input files (PE or SE).
    """
    input:
        unpack(lambda _: core.get_fq_input(config['input_type']))
    output:
        IO_FASTQ = 'fq_dict.io'
    params:
        input_type = config['input_type']
    run:
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        ReportPipeline.construct_fq_dict(input, params.input_type, Path(output.IO_FASTQ))

rule core_link_fasta_scrubbing_input:
    """
    Creates the FASTA input for the human read scrubbing step.
    """
    params:
        fasta_in = config.get('input', {}).get('fasta')
    output:
        FASTA = human_read_scrubbing.INPUT_FASTA.format(input_format='fasta')
    run:
         from camel.app.io.tooliofile import ToolIOFile
         path_fasta_in = Path(params.fasta_in[0]['path'])
         snakemakeutils.dump_object([ToolIOFile(path_fasta_in)], Path(output.FASTA))

rule core_link_vcf_input:
    """
    Creates the VCF input for the variant filtering.
    """
    params:
        vcf_in = config.get('input', {}).get('vcf_unfiltered')
    output:
        VCF = 'input/vcf.io'
    run:
         from camel.app.io.tooliofile import ToolIOFile
         path_vcf_in = Path(params.vcf_in[0]['path'])
         snakemakeutils.dump_object([ToolIOFile(path_vcf_in)], Path(output.VCF))

rule core_select_fasta:
    """
    This rules links the output of the assembly workflow to the other workflows. 
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = gene_detection.INPUT_FASTA
    shell:
        "cp {input.FASTA} {output.FASTA};"

rule core_link_fasta_to_typing:
    """
    This rule core_links the output of the assembly workflows to the sequence typing workflow.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA_typing = sequence_typing.INPUT_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA_typing};
        """

rule core_link_fasta_to_amrfinder:
    """
    This rule core_links the output of the assembly workflows to the AMRFinder workflow.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = amrfinder.INPUT_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        """

rule core_link_fasta_to_mob_suite:
    """
    This rule core_links the output of the assembly workflows to the MOB-suite workflow.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = mobsuite.INPUT_FASTA
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        """

rule core_init_summary:
    """
    Initializes the summary output file.
    """
    output:
        OUT = 'summary/summary-init.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import datetime
        from camel.app.components.pipelines.reportpipeline import ReportPipeline
        from camel.app.config import config as camel_config

        analysis_date = datetime.datetime.now().strftime(camel_config.date_fmt)
        content = [
            ('pipeline_name', config['pipeline']['name']),
            ('pipeline_version', config['pipeline']['version']),
            ('input_type', config['input_type']),
            ('sample', config['sample_name']),
            ('input_files', ReportPipeline.format_input_string(config['input'])),
            ('analysis_date', analysis_date),
            ('detection_method', config['detection_method'])
        ]
        snakemakeutils.export_summary(content, Path(output.OUT), str(params.ext))

rule core_report_pickle_citations:
    """
    This rule core_creates a pickle with a report section containing the citations.
    """
    output:
        HTML = 'report/html-citations.iob'
    params:
        citation_keys = config['citations']
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.snakemake import snakemakeutils
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        section = SnakePipelineUtils.create_citations_section(
            params.citation_keys['other'], params.citation_keys['main'])
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule core_link_fasta_to_polishing:
    """
    Links the Flye output to the short read polishing workflow (for hybrid assembly).
    """
    input:
        FASTA = assembly_flye.OUTPUT_FASTA
    output:
        FASTA = polish_assembly_short.INPUT_ASSEMBLY_FASTA.format(assembly_type='flye')
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        """

rule core_link_flye_assembly_to_medaka_input:
    """
    Links the long-read assembly to the medaka polishing.
    """
    input:
        FASTA_flye = assembly_flye.OUTPUT_FASTA
    output:
        FASTA_medaka_flye = polish_assembly_long.INPUT_ASSEMBLY_FASTA.format(assembly_type='flye')
    run:
        shutil.copyfile(input.FASTA_flye, output.FASTA_medaka_flye)

rule core_select_summary_json:
    """
    Links the JSON summary output to the output file.
    A placeholder path is used if 'output_json' is either missing from the config file or set to None.
    """
    input:
        OUT = str(core.OUTPUT_SUMMARY).format(ext='json')
    output:
        JSON = config['output_json'] if config.get('output_json') is not None else Path(config['working_dir'], 'summary', 'summary.json')
    shell:
        """
        cp {input.OUT} {output.JSON}
        """

rule core_select_summary_tsv:
    """
    Links the TSV summary output to the output file.
    """
    input:
        OUT = str(core.OUTPUT_SUMMARY).format(ext='tsv')
    output:
        JSON = config['output_tabular']
    shell:
        """
        cp {input.OUT} {output.JSON}
        """
