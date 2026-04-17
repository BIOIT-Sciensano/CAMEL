from pathlib import Path

from camel.app.core.errors import PipelineExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.piping import pipeutils
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake.step import Step
from camel.scripts.broadwgs import references
from camel.scripts.broadwgs.snakefile import alignment


rule bwa_aln_to_bam:
    """
    Align the input fastq files with BWA and pipe the output (SAM) to Samtools to create a BAM file
    """
    input:
        FASTQ = Path(config['working_dir']) / 'input' / "{input_basename}.fastq.gz.io",
        FASTA_GENOME = Path(config['working_dir']) / references.FASTA_GENOME,
    output:
        BAM = Path(config['working_dir']) / "alignment" / "mapping" / "{input_basename}.aligned.bam.io",
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / 'alignment' / "mapping" / wildcards.input_basename,
        output_file = lambda wildcards: f'{wildcards.input_basename}.aligned.unsorted.bam'
    threads: config["params_smk"]["threads_bwa"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bwa"]
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        from camel.app.tools.samtools.samtoolsview import SamtoolsView

        Path(params.working_dir).mkdir(exist_ok=True)

        bwa_mem = BWAMap()
        sam_to_bam = SamtoolsView()

        snakemakeutils.add_io_input(bwa_mem,'FASTQ_PE', Path(input.FASTQ))
        snakemakeutils.add_io_input(bwa_mem,'INDEX_GENOME_PREFIX', Path(input.FASTA_GENOME))

        bwa_mem.update_parameters(
            threads = threads,
            **config['rule_params']['alignment'][rule]
        )

        sam_to_bam.update_parameters(
            output_filename = params.output_file,
            threads = threads
        )

        pipeutils.run_as_pipe([bwa_mem, sam_to_bam], params.working_dir)

        snakemakeutils.dump_io_output(sam_to_bam,"BAM", Path(output.BAM))

rule picard_add_readgroups:
    """
    Add readgroup information to the BAM files. In Broad's WDL pipeline, the input files are in the uBAM format and 
    already contain this information. ReadGroup information is required for some of the downstream steps.
    """
    input:
        BAM = rules.bwa_aln_to_bam.output.BAM
    output:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_INTERMEDIATE_BAM
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'add_readgroups',
        output_file = lambda wildcards: f'{wildcards.input_basename}.aligned_rgadded.bam',
        RG_id = lambda wildcards: f'{wildcards.input_basename}.rg'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups

        add_rg = AddOrReplaceReadGroups()
        snakemakeutils.add_io_inputs(add_rg, input)
        step = Step(rule_name=str(rule), tool=add_rg, dir_=params.working_dir)
        add_rg.update_parameters(
            output = params.output_file,
            RG_id = params.RG_id,
            RG_sample_name = config['sample'],
            **config['rule_params']['alignment'][rule],
        )
        step.run()
        snakemakeutils.dump_io_output(add_rg,"BAM", Path(output.BAM))

rule picard_mark_duplicates_sort:
    """
    Aggregate aligned+merged flowcell BAM files and mark duplicates
    We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
    to avoid having to spend time just merging BAM files.
    """
    input:
        BAM = expand(rules.picard_add_readgroups.output.BAM, input_basename = config['input_basenames'])
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'mark_duplicates' / f'{config["sample"]}.aligned.sorted.duplicates_marked.bam.io',
        metrics = Path(config['working_dir']) / alignment.OUTPUT_MARK_DUPLICATES_METRICS
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'mark_duplicates',
        qc_dir = Path(config['working_dir']) / 'qc' / 'mark_duplicates',
        output_file = f'{config["sample"]}.aligned.unsorted.duplicates_marked.bam',
        metrics_output_file = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / "duplicate_metrics.txt"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_mark_duplicates_sort"]
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates
        from camel.app.tools.picard.sortsam import SortSam

        Path(params.qc_dir).mkdir(exist_ok=True, parents=True)

        mark_duplicates = MarkDuplicates()
        sort_sam = SortSam()

        mark_duplicates.add_input_files({"BAM": [snakemakeutils.load_object(Path(path))[0] for path in input.BAM]})

        mark_duplicates.update_parameters(
            output = params.output_file,
            metrics_output = params.metrics_output_file,
            **config['rule_params']['alignment']['picard_mark_duplicates'],
        )
        sort_sam.update_parameters(
            output = "aligned.duplicate_marked.sorted.bam",
            **config['rule_params']['alignment']['picard_sort_sam'],
        )
        mark_duplicates.update_java_options(f"-mx{resources.mem_mb}M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")
        sort_sam.update_java_options(f"-mx{resources.mem_mb}M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")

        pipeutils.run_as_pipe([mark_duplicates, sort_sam], params.working_dir)

        snakemakeutils.dump_io_output(mark_duplicates,"METRICS", Path(output.metrics))
        snakemakeutils.dump_io_output(sort_sam,"BAM", Path(output.BAM))

rule picard_set_tags:
    """
    There were discrepancies in the NM tags set by BWA, so this function recalculates them
    """
    input:
        BAM = rules.picard_mark_duplicates_sort.output.BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE
    output:
        BAM = Path(config['working_dir']) / 'alignment' / 'set_tags' / 'aligned_rgadded_dup-removed_settags.bam.io'
    params:
        working_dir = Path(config['working_dir']) / 'alignment' / 'set_tags',
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.setnmmdanduqtags import SetNmMdAndUqTags

        set_tags = SetNmMdAndUqTags()
        snakemakeutils.add_io_inputs(set_tags, input)
        step = Step(rule_name=str(rule), tool=set_tags, dir_=params.working_dir)
        set_tags.update_parameters(
            output = "aligned_dupmarked_sorted_rgadded_settags.bam",
            **config['rule_params']['alignment'][rule]
        )
        step.run()
        snakemakeutils.dump_io_output(set_tags,"BAM", Path(output.BAM))

checkpoint create_intervalfiles:
    """
    Generate interval files for everything BQSR-related. For efficiency, these steps will be run over several intervals 
    in parallel and gathered afterwards. Using a checkpoint will trigger reevaluation of the DAG, allowing for an unknown 
    number of output files. 
    Note: The number of intervals (interval files) is not likely to change in the near future, but implementing it in this
    way, makes it more clean and flexible if things would change
    """
    input:
        DICT_GENOME = Path(config['working_dir']) / references.DICT_GENOME
    output:
        TXT_intervals =  directory(Path(config['working_dir']) / "alignment" / "intervals")
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "intervals",
        output_txt = Path(config['working_dir']) / "alignment" / "intervals" / "sequence_intervals"
    run:
        Path(params.working_dir).mkdir(exist_ok=True)

        dict_genome = snakemakeutils.load_object(Path(input.DICT_GENOME))[0].path
        with open(dict_genome, "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0
            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    #Tuple: (Sequence_Name (SN), Sequence_Length(SL))
                    sequence_name = line_split[1].split("SN:")[1]
                    sequence_length = int(line_split[2].split("LN:")[1])
                    sequence_tuple_list.append((sequence_name, sequence_length))
                    if sequence_length > longest_sequence:
                        longest_sequence = sequence_length

        if len(sequence_tuple_list) == 0:
            PipelineExecutionError("Sequence tuple list empty: no intervals available")

        # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
        # the last element after a :, so we add this as a sacrificial element.
        hg38_protection_tag = ":1+"
        # Initialize the interval with the first sequence
        interval = [sequence_tuple_list[0][0] + hg38_protection_tag]
        intervals = []

        # Initialize sequence length
        temp_size = sequence_tuple_list[0][1]
        for sequence_tuple in sequence_tuple_list[1:]:
            # Maximal size of an interval should be smaller than the longest contiguous sequence in the sequence tuple
            # list, i.c. chr 1
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                interval.append(sequence_tuple[0] + hg38_protection_tag)
            else:
                intervals.append(interval)
                interval = [sequence_tuple[0] + hg38_protection_tag]
                temp_size = sequence_tuple[1]
        # Add the unmapped sequences as a separate line (from Broad wdl pipeline)
        interval.append("unmapped")
        # Add the last interval (including unmapped) to the list of intervals
        intervals.append(interval)

        # Generate the interval files
        for n, interval in enumerate(intervals):
            with open(f'{params.output_txt}_{n}.list', "w") as fh:
                fh.write("\n".join(interval))
            snakemakeutils.dump_object([ToolIOFile(Path(f'{params.output_txt}_{n}.list'))], Path(f'{params.output_txt}_{n}.list.io'))

rule gatk4_baserecalibrator:
    """
    Run baserecalibrator over the intervals created by create_intervalfiles.
    """
    input:
        BAM = rules.picard_set_tags.output.BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
        VCF_KNOWN_SNPS = Path(config['working_dir']) / references.DBSNP,
        VCF_KNOWN_INDELS = Path(config['working_dir']) / references.KNOWN_INDELS,
        TXT_intervals = Path(config['working_dir']) / "alignment" / "intervals" / "sequence_intervals_{i}.list.io"
    output:
        TXT_RecalibrationTable = Path(config['working_dir']) / "alignment" / "bqsr" / "{i}_recal_data.csv.io"
    params:
        working_dir = lambda wildcards: Path(config['working_dir']) / "alignment" / "bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_recal_data.csv'
    threads: config["params_smk"]["threads_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4baserecalibrator import GATK4BaseRecalibrator

        bqsr = GATK4BaseRecalibrator()
        snakemakeutils.add_io_inputs(bqsr, input)

        step = Step(rule_name=str(rule), tool=bqsr, dir_=params.working_dir)
        bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run()
        snakemakeutils.dump_io_output(bqsr,"TXT_RecalibrationTable", Path(output.TXT_RecalibrationTable))

def aggregate_intervals_reports(wildcards):
    """
    Input function for the rule gatk4_gather_bqsr_reports. Re-evaluated upon completion of the checkpoint.
    """
    # Ensure that snakemake records the checkpoint as direct dependency of the rule gatk4_gather_bqsr_reports
    checkpoint_output = Path(checkpoints.create_intervalfiles.get(**wildcards).output[0])
    # Retrieve values of wildcard i based on all sequence_intervals_{i}.list.io files created in the checkpoint output directory
    # to expand output of gatk4_baserecalibrator
    return expand(Path(rules.gatk4_baserecalibrator.output.TXT_RecalibrationTable),
                  i = glob_wildcards(Path.joinpath(checkpoint_output, "sequence_intervals_{i}.list.io")).i)

rule gatk4_gather_bqsr_reports:
    """
    Function to gather the BQSR reports created over several intervals.
    """
    input:
        bqsr_report_interval = aggregate_intervals_reports
    output:
        bqsr_report_gathered = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports" / "recal_data.csv.io",
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_reports",
    threads: config["params_smk"]["threads_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4gatherbqsrreports import GATK4GatherBQSRReports

        Path(params.working_dir).mkdir(exist_ok=True)

        input_bqsr_reports = aggregate_intervals_reports(wildcards)

        gather_bqsr = GATK4GatherBQSRReports()
        step = Step(rule_name=str(rule), tool=gather_bqsr, dir_=params.working_dir)
        gather_bqsr.add_input_files({"TXT_intervals": [snakemakeutils.load_object(Path(path))[0] for path in input_bqsr_reports]})
        gather_bqsr.update_parameters(
            output = "recal_data.csv"
        )
        step.run()
        snakemakeutils.dump_io_output(gather_bqsr,'TXT_RecalibrationTable', Path(output.bqsr_report_gathered))

rule gatk4_apply_bqsr:
    """
    Recalibrate the base qualities of the input reads based on the recalibration table produced by gatk4_baserecalibrator
    """
    input:
        BAM = rules.gatk4_baserecalibrator.input.BAM, # Same input as baserecalibrator
        BQSR = rules.gatk4_gather_bqsr_reports.output.bqsr_report_gathered,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
        TXT_intervals = rules.gatk4_baserecalibrator.input.TXT_intervals
    output:
        BAM = Path(config['working_dir']) / "alignment" / "apply_bqsr" / "{i}_sorted.bam.io"
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "apply_bqsr",
        output_file = lambda wildcards: f'{wildcards.i}_sorted.bam'
    threads: config["params_smk"]["threads_apply_bqsr"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_bqsr"]
    run:
        from camel.app.tools.gatk4.gatk4applybqsr import GATK4ApplyBQSR

        Path(params.working_dir).mkdir(exist_ok=True)

        apply_bqsr = GATK4ApplyBQSR()
        snakemakeutils.add_io_inputs(apply_bqsr, input)
        step = Step(rule_name=str(rule), tool=apply_bqsr, dir_=params.working_dir)
        apply_bqsr.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = params.output_file
        )
        step.run_step()
        snakemakeutils.dump_io_output(apply_bqsr,"BAM", Path(output.BAM))

def aggregate_intervals_bam(wildcards):
    """
    Input function for the rule picard_gather_sorted_bam. Re-evaluated upon completion of the checkpoint.
    """
    # Ensure that snakemake records the checkpoint as direct dependency of the rule picard_gather_sorted_bam
    checkpoint_output = Path(checkpoints.create_intervalfiles.get(**wildcards).output[0])
    # Retrieve values of wildcard i based on all sequence_intervals_{i}.list.io files created in the checkpoint output directory
    # to expand on the output of gatk4_apply_bqsr
    return expand(Path(rules.gatk4_apply_bqsr.output.BAM),
                  i = glob_wildcards(Path.joinpath(checkpoint_output, "sequence_intervals_{i}.list.io")).i)

rule picard_gather_sorted_bam:
    """
    Gather the BAM files created over multiple intervals by gatk4_apply_bqsr.
    This Base-score recalibrated, gathered, sorted BAM is the final output of the alignment step.
    """
    input:
        bqsr_BAM_interval = aggregate_intervals_bam
    output:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
    params:
        working_dir = Path(config['working_dir']) / "alignment" / "gather_bqsr_sorted_bam"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.gatherbamfiles import GatherBamFiles

        Path(params.working_dir).mkdir(exist_ok=True)

        gather_bam = GatherBamFiles()
        step = Step(rule_name=str(rule), tool=gather_bam, dir_=params.working_dir)
        gather_bam.add_input_files({"BAMs": [snakemakeutils.load_object(Path(path))[0] for path in input.bqsr_BAM_interval]})
        gather_bam.update_parameters(
            **config['rule_params']['alignment'][rule],
            output = f"{config['sample']}_gathered_sorted.bam"
        )
        step.run_step()
        snakemakeutils.dump_io_output(gather_bam,'BAM', Path(output.BAM))
