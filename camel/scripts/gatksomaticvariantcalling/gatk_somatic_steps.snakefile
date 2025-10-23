import os

from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.core.snakemake.step import Step
from camel.app.snakemake.snakemakeutils import SnakemakeUtils

working_dir = config['working_dir']
TOOL_PARAM_DIR = config['TOOL_PARAM_DIR']

camel = Camel(tool_parameter_loc=TOOL_PARAM_DIR)

def prepare_bwa_alignment_input(wildcards):
    """
    Prepares input for bwa_alignment rule. 
    Acts as a static fork in workflow based on config file (execution of first base trimming or not). 
    :return: BAM iofile path
    """
    if config["run_remove_first_base"]==True:
        FASTQ = os.path.join(working_dir, "removefirstbase/fastq.io")
    else:
        FASTQ = os.path.join(working_dir, "initial_input/fastq.io")
    return FASTQ

def prepare_addreadgroups_input(wildcards):
    """
    Prepares input for addreadgroups rule. 
    Acts as a static fork in workflow based on config file (execution of markduplicates or not). 
    :return: 
    """
    if config["run_markDuplicates"]:
        BAM = os.path.join(working_dir, "markduplicates/bam.io")
    else:
        BAM = os.path.join(working_dir, "picardsortbam/sortedbam.io")
    return BAM

def prepare_printreads_input(wildcards):
    """
    Prepares input for printreads rule. 
    Acts as a static fork in workflow based on config file (execution of indel realignment or not). 
    :return: BAM iofile path
    """
    if config["run_indel_realignment"]:
        BAM = os.path.join(working_dir, "indelrealigner/bam.io")
    else:
        BAM = os.path.join(working_dir, "addreadgroups/bam.io")
    return BAM

def prepare_basequalityrecalibration_input(wildcards):
    """
    Prepares input for basequalityrecalibration rule. 
    Acts as a static fork in workflow based on config file (execution of indel realignment or not). 
    :return: BAM iofile path
    """
    if config["run_indel_realignment"]:
        BAM = os.path.join(working_dir, "indelrealigner/bam.io")
    else:
        BAM = os.path.join(working_dir, "addreadgroups/bam.io")
    return BAM

def define_pipeline_outputs(wildcards):
    """
    Defines the expected output of the pipeline for the rule "all", depending on:
    - the variant caller(s) used 
    - if using vep or oncotator to annotate or none.
    Acts as a static fork in workflow based on config file.
    :param wildcards: 
    :return: VCF iofile path
    """
    results = dict()
    results["PDF"] = os.path.join(working_dir, "analyzecovariates/pdf.io")
    results["BAM"] = os.path.join(working_dir, "printreads/bam.io")

    if "mutect1" in config["variant_caller"]:
        results["VCF_MUTECT1"] = os.path.join(working_dir, "mutect1/vcf.io"),
        results["TXT_CALL_STATS"] = os.path.join(working_dir, "mutect1/txt.io"),

    if "mutect2" in config["variant_caller"]:
        results["VCF_MUTECT2"] = os.path.join(working_dir, "mutect2/vcf.io"),
        if 'mutect2_bam_output' in config:
            results["BAM_MUTECT2"] = os.path.join(working_dir, "mutect2/bam.io"),
            results["BAI_MUTECT2"] = os.path.join(working_dir, "mutect2/bai.io"),

    if config['run_vep']:
        results["OUT_VEP"] = os.path.join(working_dir, "vep/tsv.io"),
    if config['run_oncotator']:
        results["OUT_ONCOTATOR"] = os.path.join(working_dir, "oncotator/tsv.io"),

    return results

rule all:
    """
    This rule makes sure that all other rules are executed.
    Requires the done.flag file to be present in ./output directory. This flag is created when the move_output rule finishes.
    """
    input:
        DONE = os.path.join(working_dir, "output/done.flag"),


rule move_output:
    """
    Hard-links output of pipeline to final path.
    If using as stand-alone, hard-links to 'output' directory and renames according to pipeline parameters.
    If using Galaxy, hard-links to the galaxy-defined path.
    """
    input:
        unpack(define_pipeline_outputs),

    output:
        touch(os.path.join(working_dir, "output", "done.flag")),

    run:
        # import os
        if config['from_galaxy']:
            if 'mutect1_tab_output' in config:
                tab_init_path = SnakemakeUtils.load_object(input.TXT_CALL_STATS)[0]
                try:
                    os.link(tab_init_path.path, config['mutect1_tab_output'])
                except FileExistsError:
                    os.remove(config['mutect1_tab_output'])
                    os.link(tab_init_path.path, config['mutect1_tab_output'])
            if 'mutect1_vcf_output' in config:
                vcf_init_path = SnakemakeUtils.load_object(input.VCF_MUTECT1)[0]
                try:
                    os.link(vcf_init_path.path, config['mutect1_vcf_output'])
                except FileExistsError:
                    os.remove(config['mutect1_vcf_output'])
                    os.link(vcf_init_path.path, config['mutect1_vcf_output'])
            if 'mutect2_vcf_output' in config:
                vcf_init_path = SnakemakeUtils.load_object(input.VCF_MUTECT2)[0]
                try:
                    os.link(vcf_init_path.path, config['mutect2_vcf_output'])
                except FileExistsError:
                    os.remove(config['mutect2_vcf_output'])
                    os.link(vcf_init_path.path, config['mutect2_vcf_output'])
            if 'mutect2_bam_output' in config:
                bam_init_path = SnakemakeUtils.load_object(input.BAM_MUTECT2)[0]
                try:
                    os.link(bam_init_path.path, config['mutect2_bam_output'])
                except FileExistsError:
                    os.remove(config['mutect2_bam_output'])
                    os.link(bam_init_path.path, config['mutect2_bam_output'])
            if 'bam_output' in config:
                bam_init_path = SnakemakeUtils.load_object(input.BAM)[0]
                try:
                    os.link(bam_init_path.path, config['bam_output'])
                except FileExistsError:
                    os.remove(config['bam_output'])
                    os.link(bam_init_path.path, config['bam_output'])
            if config['run_vep']:
                vep_out_init_path = SnakemakeUtils.load_object(input.OUT_VEP)[0]
                try:
                    os.link(vep_out_init_path.path, config['out_vep'])
                except FileExistsError:
                    os.remove(config['out_vep'])
                    os.link(vep_out_init_path.path, config['out_vep'])
            if config['run_oncotator']:
                tsv_oncotator_init_path = SnakemakeUtils.load_object(input.OUT_ONCOTATOR)[0]
                try:
                    os.link(tsv_oncotator_init_path.path, config['out_oncotator'])
                except FileExistsError:
                    os.remove(config['out_oncotator'])
                    os.link(tsv_oncotator_init_path.path, config['out_oncotator'])


        else:
            if 'mutect1_tab_output' in config:
                tab_init_path = SnakemakeUtils.load_object(input.TXT_CALL_STATS)[0]
                os.link(tab_init_path.path, os.path.join(working_dir, "output/", config['mutect1_tab_output']))
            if 'mutect1_vcf_output' in config:
                vcf_init_path = SnakemakeUtils.load_object(input.VCF_MUTECT1)[0]
                os.link(vcf_init_path.path, os.path.join(working_dir, "output/", config['mutect1_vcf_output']))
            if 'mutect2_vcf_output' in config:
                vcf_init_path = SnakemakeUtils.load_object(input.VCF_MUTECT2)[0]
                os.link(vcf_init_path.path, os.path.join(working_dir, "output/", config['mutect2_vcf_output']))
            if 'mutect2_bam_output' in config:
                bam_init_path = SnakemakeUtils.load_object(input.BAM_MUTECT2)[0]
                os.link(bam_init_path.path, os.path.join(working_dir, "output/", config['mutect2_bam_output']))
                bai_init_path = SnakemakeUtils.load_object(input.BAI_MUTECT2)[0]
                if os.path.splitext(config['mutect2_bam_output'])[1] == ".bam":
                    bai_output = os.path.splitext(config['mutect2_bam_output'])[0]+".bai"
                else:
                    bai_output = config['mutect2_bam_output']+".bai"
                os.link(bai_init_path.path, os.path.join(working_dir, "output/", bai_output))
            if 'covar_output' in config:
                pdf_init_path = SnakemakeUtils.load_object(input.PDF)[0]
                os.link(pdf_init_path.path, os.path.join(working_dir, "output/", config['covar_output']))
            if 'bam_output' in config:
                bam_init_path = SnakemakeUtils.load_object(input.BAM)[0]
                os.link(bam_init_path.path, os.path.join(working_dir, "output/", config['bam_output']))
            if config['run_vep']:
                out_vep_init_path = SnakemakeUtils.load_object(input.OUT_VEP)[0]
                os.link(out_vep_init_path.path, os.path.join(working_dir, "output/", config['out_vep']))
            if config['run_oncotator']:
                tsv_oncotator_init_path = SnakemakeUtils.load_object(input.OUT_ONCOTATOR)[0]
                os.link(tsv_oncotator_init_path.path, os.path.join(working_dir, "output/", config['out_oncotator']))


rule prepare_initial_input:
    """ Prepare input for the pipeline: generates io file from PE or SE fastq file(s) """
    input:
        FASTQ=config['fastq'],
    output:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
    run:
        SnakemakeUtils.pickle_snake_input(input, output)



rule prepare_references_io:
    """
    Prepare reference genome IO files for snakemake to use: generate io files for
    - the reference genome fasta (value and file objects)
    - the snp vcf file
    - the indel vcf file.
    Requires reference names from db_loc to be present in the config file.
    """
    output:
        FASTA_GENOME=os.path.join(working_dir, "initial_input/fasta_reference_human_value.io"),
        FASTA_GENOME_FILE=os.path.join(working_dir, "initial_input/fasta_reference_human.io"), # for tools that require ToolIOFile instead of ToolIOValue
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    run:
        from camel.app.io.tooliodb import ToolIODb
        fasta_genome = config['fasta_ref'],
        vcf_known_snps = config['vcf_known_snps'],
        vcf_known_indels = config['vcf_known_indels'],
        io_fasta_genome = [ToolIOValue(ToolIODb(fasta_genome).path)]
        io_file_fasta_genome = [ToolIOFile(ToolIODb(fasta_genome).path)]
        io_vcf_known_snps = [ToolIOFile(ToolIODb(vcf_known_snps).path)]
        io_vcf_known_indels = [ToolIOFile(ToolIODb(vcf_known_indels).path)]
        SnakemakeUtils.dump_object(io_fasta_genome, output.FASTA_GENOME)
        SnakemakeUtils.dump_object(io_file_fasta_genome, output.FASTA_GENOME_FILE)
        SnakemakeUtils.dump_object(io_vcf_known_snps, output.VCF_KNOWN_SNPS)
        SnakemakeUtils.dump_object(io_vcf_known_indels, output.VCF_KNOWN_INDELS)


rule removefirstbases:
    """
    Remove first base of reads to avoid technical artifacts in amplicon sequencing.
    """
    input:
        FASTQ=os.path.join(working_dir, "initial_input/fastq.io"),
    output:
        FASTQ=os.path.join(working_dir, "removefirstbase/fastq.io")
    threads: 8
    params:
        working_dir=os.path.join(working_dir, "removefirstbase"),
    run:
        from camel.app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimm = Trimmomatic(camel)
        if config['PE']:
            SnakemakeUtils.add_pickle_input(trimm, 'FASTQ_PE', input.FASTQ)
        if config['SE']:
            SnakemakeUtils.add_pickle_input(trimm, 'FASTQ_SE', input.FASTQ)
        trimm.update_parameters(threads=threads)
        trimm.update_parameters(head_crop=config['bases_to_rm_from_reads'])
        step = Step(rule_name=str(rule), tool=trimm, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(trimm, "FASTQ_PE", output.FASTQ)



rule bwa_alignment:
    """
    Reads alignment using bwa mem.
    """
    input:
        FASTQ=prepare_bwa_alignment_input,
        FASTA_GENOME=os.path.join(working_dir, "initial_input/fasta_reference_human_value.io"),
    output:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io")
    threads: 16
    params:
        working_dir=os.path.join(working_dir, "bwa_alignment"),
    run:
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_mem = BWAMap(camel)
        if config['PE']:
            SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_PE', input.FASTQ)
        if config['SE']:
            SnakemakeUtils.add_pickle_input(bwa_mem, 'FASTQ_SE', input.FASTQ)
        SnakemakeUtils.add_pickle_input(bwa_mem, 'INDEX_GENOME_PREFIX', input.FASTA_GENOME)
        step = Step(rule_name=str(rule), tool=bwa_mem, params.working_dir, config)
        bwa_mem.update_parameters(threads=threads)
        step.run()
        SnakemakeUtils.dump_tool_output(bwa_mem, "SAM", output.SAM)


rule samtobam:
    """
    Sam-to-bam conversion (samtools).
    """
    input:
        SAM=os.path.join(working_dir, "bwa_alignment/sam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "bwa_alignment")
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        smv = SamtoolsView(camel)
        SnakemakeUtils.add_pickle_input(smv, "SAM", input.SAM)
        step = Step(rule_name=str(rule), tool=smv, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(smv, "BAM", output.BAM)


rule sortbam:
    """
    pre-indexing Bam sorting (samtools). 
    """
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/bam.io"),
    output:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    params:
        working_dir = os.path.join(working_dir, "bwa_alignment")
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        sms = SamtoolsSort(camel)
        SnakemakeUtils.add_pickle_input(sms,"BAM",input.BAM)
        step = Step(rule_name=str(rule), tool=sms, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(sms, "BAM", output.BAM)


rule indexbam:
    """
    Bam indexing (samtools).
    """
    input:
        BAM=os.path.join(working_dir, "bwa_alignment/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    params:
        working_dir = os.path.join(working_dir, "samtools_index")
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        smi = SamtoolsIndex(camel)
        SnakemakeUtils.add_pickle_input(smi,"BAM",input.BAM)
        step = Step(rule_name=str(rule), tool=smi, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(smi, "BAM", output.BAM)


rule picardsortbam:
    """
    Bam sorting (picard).
    """
    input:
        BAM=os.path.join(working_dir, "samtools_index/bam_after_index.io"),
    output:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    params:
        working_dir = os.path.join(working_dir, "picardsortbam")
    run:
        from camel.app.tools.picard.sortsam import SortSam
        pss = SortSam(camel)
        SnakemakeUtils.add_pickle_input(pss,"BAM",input.BAM)
        step = Step(rule_name=str(rule), tool=pss, params.working_dir, config)
        pss.update_parameters(sort_order='coordinate')
        step.run()
        SnakemakeUtils.dump_tool_output(pss, "BAM", output.BAM)


rule markduplicates:
    """
    Optional reads deduplication (picard).
    """
    input:
        BAM=os.path.join(working_dir, "picardsortbam/sortedbam.io"),
    output:
        BAM=os.path.join(working_dir, "markduplicates/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "markduplicates")
    run:
        from camel.app.tools.picard.markduplicates import MarkDuplicates
        pmd = MarkDuplicates(camel)
        SnakemakeUtils.add_pickle_input(pmd, "BAM", input.BAM)
        step = Step(rule_name=str(rule), tool=pmd, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(pmd, "BAM", output.BAM)


rule addreadgroups:
    """
    Read group adding (picard)
    """
    input:
        BAM=prepare_addreadgroups_input,
    output:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "addreadgroups")
    run:
        from camel.app.tools.picard.addorreplacereadgroups import AddOrReplaceReadGroups
        parg = AddOrReplaceReadGroups(camel)
        SnakemakeUtils.add_pickle_input(parg, "BAM", input.BAM)
        step = Step(rule_name=str(rule), tool=parg, params.working_dir, config)
        parg.update_parameters(create_index='true')
        step.run()
        SnakemakeUtils.dump_tool_output(parg, "BAM", output.BAM)


rule bamtobed:
    """
    Bam-to-bed conversion for interval generation for optimisation (bedtools).
    """
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
    output:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    params:
        working_dir = os.path.join(working_dir, "bamtobed")
    run:
        from camel.app.tools.bedtools.bedtoolsbamtobed import BedtoolsBamToBed
        btb = BedtoolsBamToBed(camel)
        SnakemakeUtils.add_pickle_input(btb, "BAM", input.BAM)
        step = Step(rule_name=str(rule), tool=btb, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(btb,"BED",output.BED)


rule generate_intervals:
    """
    Interval generation for optimisation (bedtools).
    """
    input:
        BED=os.path.join(working_dir, "bamtobed/bed.io"),
    output:
        BED=os.path.join(working_dir, "generate_intervals/bed.io"),
    params:
        working_dir = os.path.join(working_dir, "generate_intervals")
    run:
        from camel.app.tools.bedtools.bedtoolsmerge import BedtoolsMerge
        btm = BedtoolsMerge(camel)
        SnakemakeUtils.add_pickle_input(btm, "BED", input.BED)
        step = Step(rule_name=str(rule), tool=btm, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(btm, "BED", output.BED)


rule realignertargetcreator:
    """
    Indel realigner intervals creation (GATK).
    """
    input:
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        TXT_intervals=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),

    output:
        INTERVALS=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
    params:
        working_dir = os.path.join(working_dir, "realignertargetcreator")
    run:
        from camel.app.tools.gatk.gatkrealignertargetcreator import GATKRealignerTargetCreator
        grtc = GATKRealignerTargetCreator(camel)
        SnakemakeUtils.add_pickle_inputs(grtc, input)
        step = Step(rule_name=str(rule), tool=grtc, params.working_dir, config)
        step.run()
        SnakemakeUtils.dump_tool_output(grtc, "TXT_realign_intervals", output.INTERVALS)


rule indelrealigner:
    """
    Indel realignment (GATK).
    """
    input:
        TXT_realign_intervals=os.path.join(working_dir, "realignertargetcreator/intervals.io"),
        BAM=os.path.join(working_dir, "addreadgroups/bam.io"),
        TXT_intervals=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
    output:
        BAM=os.path.join(working_dir, "indelrealigner/bam.io"),
    params:
        working_dir = os.path.join(working_dir, "indelrealigner"),
    run:
        from camel.app.tools.gatk.gatkindelrealigner import GATKIndelRealigner
        gir=GATKIndelRealigner(camel)
        SnakemakeUtils.add_pickle_inputs(gir, input)
        step = Step(rule_name=str(rule), tool=gir, params.working_dir, config)
        step.run_step()
        SnakemakeUtils.dump_tool_output(gir,"BAM",output.BAM)


rule basequalityrecalibration:
    """
    Base quality recalibration (GATK).
    """
    input:
        BAM=prepare_basequalityrecalibration_input,
        TXT_intervals=os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    output:
        TXT=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "basequalityrecalibration"),
    run:
        from camel.app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        bqsr=GATKBaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_inputs(bqsr, input)
        step = Step(rule_name=str(rule), tool=bqsr, params.working_dir, config)
        bqsr.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr,"TXT_RecalibrationTable",output.TXT)


rule printreads:
    """
    Post-recalibration reads printing (GATK). 
    """
    input:
        BAM=prepare_printreads_input,
        BQSR=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        TXT_intervals = os.path.join(working_dir, "generate_intervals/bed.io"),
    output:
        BAM=os.path.join(working_dir, "printreads/bam.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "printreads"),
    run:
        from camel.app.tools.gatk.gatkprintreads import GATKPrintReads
        gpr=GATKPrintReads(camel)
        SnakemakeUtils.add_pickle_inputs(gpr, input)
        step = Step(rule_name=rule, tool=gpr, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        gpr.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_output(gpr,"BAM",output.BAM)


rule basequalityrecalibration2:
    """
    Second base quality recalibration for quality control (GATK).
    """
    input:
        BAM = os.path.join(working_dir, "printreads/bam.io"),
        TXT_intervals = os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_KNOWN_SNPS=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
        VCF_KNOWN_INDELS=os.path.join(working_dir, "initial_input/vcf_known_indels.io"),
    output:
        TXT = os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
    threads: 5
    params:
        working_dir = os.path.join(working_dir, "basequalityrecalibration2"),
    run:
        from camel.app.tools.gatk.gatkbaserecalibrator import GATKBaseRecalibrator
        bqsr = GATKBaseRecalibrator(camel)
        SnakemakeUtils.add_pickle_inputs(bqsr, input)
        step = Step(rule_name=str(rule), tool=bqsr, params.working_dir, config)
        bqsr.update_parameters(threads=threads)
        step.run_step()
        SnakemakeUtils.dump_tool_output(bqsr, "TXT_RecalibrationTable", output.TXT)


rule analyzecovariates:
    """
    Covariates analysis and BQSR report generation (GATK) 
    """
    input:
        TXT_TABLE_BEFORE=os.path.join(working_dir, "basequalityrecalibration/txt.io"),
        TXT_TABLE_AFTER=os.path.join(working_dir, "basequalityrecalibration2/txt.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
    output:
        PDF=os.path.join(working_dir, "analyzecovariates/pdf.io"),
    params:
        working_dir = os.path.join(working_dir, "analyzecovariates")
    run:
        from camel.app.tools.gatk.gatkanalyzecovariates import GATKAnalyzeCovariates
        gac = GATKAnalyzeCovariates(camel)
        SnakemakeUtils.add_pickle_inputs(gac, input)
        step = Step(rule_name=rule, tool=gac, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        if 'covar_output' in config:
            gac.update_parameters(pdf_output = config['covar_output'])
        step.run_step()
        SnakemakeUtils.dump_tool_output(gac, "PDF", output.PDF)


rule mutect1:
    """
    Variant calling (MuTect1).
    """
    input:
        BAM_TUMOR=os.path.join(working_dir, "printreads/bam.io"),
        TXT_intervals = os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_DBSNP=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
    output:
        TXT_CALL_STATS=os.path.join(working_dir, "mutect1/txt.io"),
        VCF=os.path.join(working_dir, "mutect1/vcf.io"),
    params:
        working_dir = os.path.join(working_dir, "mutect1"),
    run:
        from camel.app.tools.mutect.mutect1 import Mutect1
        mut=Mutect1(camel)
        SnakemakeUtils.add_pickle_inputs(mut, input)
        step = Step(rule_name=rule, tool=mut, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        if 'MuTect1_downsampling_target' in config:
            mut.update_parameters(downsampling_coverage_target=config['MuTect1_downsampling_target'])
        if 'MuTect1_downsampling_type' in config:
            mut.update_parameters(downsampling_type=config['MuTect1_downsampling_type'])
        if 'gap_events_threshold' in config:
            mut.update_parameters(gap_events_threshold=config['gap_events_threshold'])
        if 'strand_artifact_lod' in  config:
            mut.update_parameters(strand_artifact_lod=config['strand_artifact_lod'])
        step.run_step()
        SnakemakeUtils.dump_tool_outputs(mut, output)


rule mutect2:
    """
    Variant calling (MuTect2).
    """
    input:
        BAM_TUMOR=os.path.join(working_dir, "printreads/bam.io"),
        TXT_intervals = os.path.join(working_dir, "generate_intervals/bed.io"),
        FASTA_REF=os.path.join(working_dir, "initial_input/fasta_reference_human.io"),
        VCF_DBSNP=os.path.join(working_dir, "initial_input/vcf_known_snps.io"),
    output:
        VCF=os.path.join(working_dir, "mutect2/vcf.io"),
        BAM=os.path.join(working_dir, "mutect2/bam.io"),
        BAI=os.path.join(working_dir, "mutect2/bai.io"),
    params:
        working_dir = os.path.join(working_dir, "mutect2"),
    run:
        from camel.app.tools.gatk.gatkmutect2 import GATKMuTect2
        mut2=GATKMuTect2(camel)
        SnakemakeUtils.add_pickle_inputs(mut2, input)
        step = Step(rule_name=rule, tool=mut2, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        if 'mutect_nct' in config:
            mut2.update_parameters(threads=config['mutect_nct'])
        if 'mutect2_bam_output' in config:
            mut2.update_parameters(output_bam=True)
        if 'MuTect2_downsampling_target' in config:
            mut2.update_parameters(downsampling_coverage_target=config['MuTect2_downsampling_target'])
        # Testing
        if 'MuTect2_downsampling_type' in config:
            mut2.update_parameters(downsampling_type=config['MuTect2_downsampling_type'])
        if 'MuTect2_force_active' in config and config["MuTect2_force_active"]:
            mut2.update_parameters(force_active=True)
        if 'MuTect2_output_active_region_igv' in config and config['MuTect2_output_active_region_igv']:
            mut2.update_parameters(output_active_region_igv=True)
        if 'MuTect2_active_region_igv_file' in config:
            mut2.update_parameters(active_region_igv_file=config['MuTect2_active_region_igv_file'])
        if 'MuTect2_output_mode' in config:
            mut2.update_parameters(output_mode=config['MuTect2_output_mode'])
        if 'MuTect2_disable_optimizations' in config:
            mut2.update_parameters(disable_optimizations=True)
        step.run_step()
        # set output: bam optional, vcf always generated
        SnakemakeUtils.dump_tool_output(mut2, "VCF", output.VCF)
        if 'mutect2_bam_output' in config:
            SnakemakeUtils.dump_tool_output(mut2, "BAM", output.BAM)
            SnakemakeUtils.dump_tool_output(mut2, "BAI", output.BAI)
        else:
            SnakemakeUtils.dump_object("placeholder for empty output.", output.BAM)
            SnakemakeUtils.dump_object("placeholder for empty output.", output.BAI)

rule vep:
    """
    Variant annotation with Vep (ensembl).
    """
    input:
        # TODO: add fork if mutect2
        VCF=os.path.join(working_dir, "mutect1/vcf.io"),
    output:
        OUT=os.path.join(working_dir, "vep/tsv.io"),
    params:
        working_dir = os.path.join(working_dir, "vep"),
    run:
        from camel.app.tools.vep.vep import Vep
        from camel.app.io.tooliodb import ToolIODb
        DB_PATH = ToolIODb("Vep_homo_sapiens_93_GRCh37").path
        run_vep =Vep(camel)
        SnakemakeUtils.add_pickle_inputs(run_vep, input)
        run_vep.update_parameters(db_loc=DB_PATH)
        step = Step(rule_name=rule, tool=run_vep, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        step.run_step()
        SnakemakeUtils.dump_tool_output(run_vep, "OUT", output.OUT)


rule oncotator:
    """
    Variant annotation with oncotator (broad).
    """
    input:
        # TODO: add fork if mutect2
        IN_FILE=os.path.join(working_dir, "mutect1/vcf.io"),
    output:
        TSV=os.path.join(working_dir, "oncotator/tsv.io"),
    params:
        working_dir = os.path.join(working_dir, "oncotator"),
    run:
        from camel.app.tools.oncotator.oncotator import Oncotator
        from camel.app.io.tooliodb import ToolIODb

        DB_PATH = ToolIODb("oncotator_v1_ds_April052016").path
        CANONICAL_VARS_PATH = ToolIODb("tx_exact_uniprot_matchesAKT1_CRLF2_FGFR1").path
        octr =Oncotator(camel)
        SnakemakeUtils.add_pickle_inputs(octr, input)
        if "oncotator_tx_mode" in config:
            octr.update_parameters(tx_mode=config["oncotator_tx_mode"])
        octr.update_parameters(db_dir=DB_PATH, canonical_tx_file=CANONICAL_VARS_PATH, input_format="VCF")
        step = Step(rule_name=rule, tool=octr, camel=camel, folder=params.working_dir, config=config, pipeline_output=True)
        step.run_step()
        SnakemakeUtils.dump_tool_output(octr, "OUT_FILE", output.TSV)