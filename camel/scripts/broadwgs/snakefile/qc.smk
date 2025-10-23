from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.broadwgs import references
from camel.scripts.broadwgs.snakefile import alignment, bam_to_cram, variant_calling


rule picard_quality_yield:
    """
    Collect metrics about reads that pass quality thresholds and Illumina-specific filters.
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_INTERMEDIATE_BAM,
    output:
        TXT = Path(config['working_dir']) / "qc" / "quality_yield" / "{input_basename}.unmapped.quality_yield_metrics.io"
    params:
        working_dir = Path(config['working_dir']) / "qc" / "quality_yield",
        output_file = lambda wildcards: f"{wildcards.input_basename}.unmapped.quality_yield_metrics"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectqualityyieldmetrics import CollectQualityYieldMetrics

        quality_yield = CollectQualityYieldMetrics()
        snakemakeutils.add_pickle_inputs(quality_yield, input)
        step = Step(rule_name=str(rule), tool=quality_yield, dir_=params.working_dir)
        quality_yield.update_parameters(**config['rule_params']['qc'][rule], output = params.output_file)
        step.run()
        snakemakeutils.dump_tool_outputs(quality_yield, output)

rule picard_unsorted_RG_quality:
    """
    Collect multiple QC metrics for aligned, unsorted read-groups.
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_INTERMEDIATE_BAM,
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "unsorted_RG_quality" / "{input_basename}.unsorted_readgroup"),
                ".base_distribution_by_cycle.pdf",
                ".base_distribution_by_cycle_metrics",
                ".quality_by_cycle.pdf",
                ".quality_by_cycle_metrics",
                ".insert_size_histogram.pdf",
                ".insert_size_metrics",
                ".quality_distribution.pdf",
                ".quality_distribution_metrics"
                )
    params:
        working_dir = Path(config['working_dir']) / "qc" / "unsorted_RG_quality",
        output_prefix = lambda wildcards: f"{wildcards.input_basename}.unsorted_readgroup"
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        unsorted_RG_quality = CollectMultipleMetrics()
        snakemakeutils.add_pickle_inputs(unsorted_RG_quality, input)
        step = Step(rule_name=str(rule), tool=unsorted_RG_quality, dir_=params.working_dir)
        unsorted_RG_quality.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run()

rule picard_RG_quality:
    """
    Collect multiple QC metrics for aligned, sorted read-groups.
    """
    input:
        BAM = Path(config["working_dir"]) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "RG_quality" / f'{config["sample"]}.readgroup'),
        ".alignment_summary_metrics",
        ".gc_bias.detail_metrics",
        ".gc_bias.pdf",
        ".gc_bias.summary_metrics")
    params:
        working_dir = Path(config['working_dir']) / "qc" / "RG_quality",
        output_prefix = f'{config["sample"]}.readgroup'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        RG_quality = CollectMultipleMetrics()
        snakemakeutils.add_pickle_inputs(RG_quality, input)
        step = Step(rule_name=str(rule), tool=RG_quality, dir_=params.working_dir)
        RG_quality.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run()

# TODO maybe implement this later
# rule check_prevalidation:
#     input:
#         duplicate_metrics = Path(config['working_dir']) / alignment.OUTPUT_MARK_DUPLICATES_METRICS,
#         chimerism_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f"{config['sample']}.agg.alignment_summary_metrics",
#     output:
#         duplication_rate = "duplication_value.txt",
#         chimerism_rate = "chimerism_value.txt"
#     params:
#         working_dir = Path(config['working_dir']) / "qc" / "check_prevalidation",
#         max_duplication_in_reasonable_sample = 0.30,
#         max_chimerism_in_reasonable_sample = 0.15,
#         #Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
#     run:
#         import subprocess
#
#         subprocess.run(f"grep -A 1 PERCENT_DUPLICATION {input.duplication_metrics} > {params.working_dir}.duplication.csv")
#         subprocess.run(f"grep -A 3 PCT_CHIMERAS {input.chimerism_metrics} | grep -v OF_PAIR > {params.working_dir}.chimerism.csv")
#
#         import csv
#         with open(f"{params.working_dir}.duplication.csv") as dupfile:
#             reader = csv.DictReader(dupfile, delimiter='\t')
#             for row in reader:
#                 with open(f"{params.working_dir}/{output.duplication_rate}","w") as file:
#                     file.write(row['PERCENT_DUPLICATION'])
#                     file.close()
#
#         with open(f"{params.working_dir}.chimerism.csv") as chimfile:
#             reader = csv.DictReader(chimfile, delimiter='\t')
#             for row in reader:
#                 with open(f"{params.working_dir}/{output.chimerism_rate}","w") as file:
#                     file.write(row['PCT_CHIMERAS'])
#                     file.close()


rule picard_aggregation_metrics:
    """
    Collect aggregation metrics
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
    output:
        multiext(str(Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg'),
        ".alignment_summary_metrics",
        ".bait_bias_detail_metrics",
        ".bait_bias_summary_metrics",
        ".gc_bias.detail_metrics",
        ".gc_bias.pdf",
        ".gc_bias.summary_metrics",
        ".insert_size_histogram.pdf",
        ".insert_size_metrics",
        ".pre_adapter_detail_metrics",
        ".pre_adapter_summary_metrics",
        ".quality_distribution.pdf",
        ".quality_distribution_metrics",
        ".error_summary_metrics")
    params:
          working_dir = Path(config['working_dir']) / "qc" / "aggregation_metrics",
          output_prefix = f'{config["sample"]}.agg',
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectmultiplemetrics import CollectMultipleMetrics

        agg_metrics = CollectMultipleMetrics()
        snakemakeutils.add_pickle_inputs(agg_metrics, input)
        step = Step(rule_name=str(rule), tool=agg_metrics, dir_=params.working_dir)
        agg_metrics.update_parameters(
            output_prefix = params.output_prefix,
            **config['rule_params']['qc'][rule]
        )
        step.run()

rule picard_RG_checksum:
    """
    Create a hash code based on the read groups (RG).
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
    output:
        TXT_checksum = Path(config['working_dir']) / "qc" / "RG_checksum" / f'{config["sample"]}.bam.read_group_md5',
    params:
        working_dir = Path(config['working_dir']) / "qc" / "RG_checksum",
        output_file = f'{config["sample"]}.bam.read_group_md5'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.calculatereadgroupchecksum import CalculateReadGroupChecksum

        get_checksum = CalculateReadGroupChecksum()
        snakemakeutils.add_pickle_input(get_checksum, "BAM", Path(input.BAM))
        step = Step(rule_name=str(rule), tool=get_checksum, dir_=params.working_dir)
        get_checksum.update_parameters(output = params.output_file)
        step.run()


rule picard_wgs_metrics:
    """
    Collect metrics about coverage and performance of WGS experiments.
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
        COVERAGE_INTERVALS = Path(config['working_dir']) / references.COVERAGE_INTERVALS
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.wgs.metrics.txt'
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        output_file = f'{config["sample"]}.wgs.metrics.txt'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectwgsmetrics import CollectWgsMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        wgs_metrics = CollectWgsMetrics()
        snakemakeutils.add_pickle_inputs(wgs_metrics, input)
        step = Step(rule_name=str(rule), tool=wgs_metrics, dir_=params.working_dir)
        wgs_metrics.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run()

rule picard_raw_wgs_metrics:
    """
    Collect whole genome sequencing-related metrics
    """
    input:
        BAM = Path(config['working_dir']) / alignment.OUTPUT_ALIGNMENT_BAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE,
        COVERAGE_INTERVALS = Path(config['working_dir']) / references.COVERAGE_INTERVALS
    output:
        TXT_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.raw.wgs.metrics.txt'
    params:
        working_dir = Path(config['working_dir']) / "qc" / "wgs_metrics",
        output_file = f'{config["sample"]}.raw.wgs.metrics.txt'
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectrawwgsmetrics import CollectRawWgsMetrics

        raw_wgs_metrics = CollectRawWgsMetrics()
        snakemakeutils.add_pickle_inputs(raw_wgs_metrics, input)
        step = Step(rule_name=str(rule), tool=raw_wgs_metrics, dir_=params.working_dir)
        raw_wgs_metrics.update_parameters(
            output = params.output_file,
            **config['rule_params']['qc'][rule]
        )
        step.run()


rule picard_validate_cram:
    """
    Validate the format of CRAM file
    """
    input:
        CRAM = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM,
        FASTA_REF = Path(config['working_dir']) / references.FASTA_GENOME_FILE
    output:
        TXT_metrics = Path(config['working_dir']) / bam_to_cram.OUTPUT_BAMTOCRAM_CRAM_metrics,
    params:
        working_dir = Path(config['working_dir']) / "qc" / "bamtocram"
    threads: config["params_smk"]["threads_cram"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_cram"]
    run:
        from camel.app.tools.picard.validatesamfile import ValidateSamFile

        Path(params.working_dir).mkdir(exist_ok=True)

        val_cram = ValidateSamFile()
        snakemakeutils.add_pickle_input(val_cram, "BAM", Path(input.CRAM))
        snakemakeutils.add_pickle_input(val_cram, "FASTA_REF", Path(input.FASTA_REF))
        step = Step(rule_name=str(rule), tool=val_cram, dir_=params.working_dir)
        val_cram.update_parameters(
            **config['rule_params']['bam_to_cram'][rule]
        )
        val_cram.update_java_options("-mx100G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false")
        step.run()
        snakemakeutils.dump_tool_output(val_cram, 'TXT_report', Path(output.TXT_metrics))

rule picard_variant_calling_metrics:
    """
    Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file.
    """
    input:
        VCF = Path(config['working_dir']) / variant_calling.OUTPUT_gVCF,
        VCF_dbsnp = Path(config['working_dir']) / references.DBSNP,
        DICT_GENOME = Path(config['working_dir']) / references.DICT_GENOME,
        EVALUATION_INTERVALS = Path(config['working_dir']) / references.EVALUATION_INTERVALS,
    output:
        TXT_report = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f'{config["sample"]}.variant_calling_metrics.io',
    params:
        working_dir = Path(config['working_dir']) / "qc" / "variant_calling_metrics",
    threads: config["params_smk"]["threads_picard"]
    resources:
        mem_mb=config["params_smk"]["memory_mb_picard"]
    run:
        from camel.app.tools.picard.collectvariantcallingmetrics import CollectVariantCallingMetrics

        Path(params.working_dir).mkdir(exist_ok=True)

        vcf_metrics = CollectVariantCallingMetrics()
        snakemakeutils.add_pickle_inputs(vcf_metrics, input)
        step = Step(rule_name=str(rule), tool=vcf_metrics, dir_=params.working_dir)
        vcf_metrics.update_parameters(
            output_prefix = config["sample"],
            **config['rule_params']['qc'][rule]
        )
        step.run()
        snakemakeutils.dump_tool_outputs(vcf_metrics, output)
