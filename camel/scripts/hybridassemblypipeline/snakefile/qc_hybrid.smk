from pathlib import Path

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.core.snakemake.step import Step
from camel.scripts.hybridassemblypipeline.snakefile import qc_hybrid


rule qc_hybrid_samtools_index:
    """
    Creates a samtools index for the assembly.
    """
    input:
        FASTA = lambda wildcards: qc_hybrid.consensus_by_tool[wildcards.name]
    output:
        FASTA = 'qc_hybrid/{name}/fasta-index.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}'
    run:
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools = SamtoolsFastaIndex()
        snakemakeutils.add_pickle_inputs(samtools, input)
        step = Step(str(rule), samtools, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools, output)

rule qc_hybrid_quast:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        TSV = 'qc_hybrid/{name}/quast/work/report.tsv',
        INFORMS = 'qc_hybrid/{name}/quast/tool/commands.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/quast'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast()
        snakemakeutils.add_pickle_inputs(quast, input)
        step = Step(str(rule), quast, dir_=Path(str(params.dir_)))
        step.run()
        quast.informs['_tag'] = wildcards.name
        snakemakeutils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_hybrid_parse_quast_output:
    """
    Parses the quast output into a pickle.
    """
    input:
        TSV = rules.qc_hybrid_quast.output.TSV
    output:
        INFORMS = 'qc_hybrid/{name}/quast/informs.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/quast'
    run:
        from camel.app.tools.quast.quastinformextractor import QuastInformExtractor
        quast_inform_extractor = QuastInformExtractor()
        quast_inform_extractor.add_input_files({'TSV': [ToolIOFile(Path(input.TSV))]})
        step = Step(str(rule), quast_inform_extractor, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(quast_inform_extractor, output)

rule qc_hybrid_quast_all_assemblies:
    """
    Generates assembly statistics using QUAST.
    """
    input:
        FASTA = [f'qc_hybrid/{name}/fasta-index.io' for name in config['assembly_steps']]
    output:
        TSV = 'qc_hybrid/quast_combined/work/report.tsv',
        HTML = 'qc_hybrid/quast_combined/work/report.html',
        INFORMS = 'qc_hybrid/quast_combined/tool/commands.io'
    params:
        dir_ = 'qc_hybrid/quast_combined'
    run:
        from camel.app.tools.quast.quast import Quast
        quast = Quast()
        quast.add_input_files({'FASTA': [snakemakeutils.load_object(Path(str(i)))[0] for i in input.FASTA]})
        step = Step(str(rule), quast, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_object(quast.informs, Path(output.INFORMS))

rule qc_hybrid_bwa_index:
    """
    Creates a bwa index for the assembly.
    """
    input:
        FASTA_REF = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        INDEX_GENOME_PREFIX = 'qc_hybrid/{name}/genome_prefix.iob'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}'
    run:
        from camel.app.tools.bwa.bwaindex import BWAIndex
        bwa = BWAIndex()
        snakemakeutils.add_pickle_inputs(bwa, input)
        step = Step(str(rule), bwa, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bwa, output)

rule qc_hybrid_read_mapping_illumina:
    """
    Maps the short reads against the assembly.
    """
    input:
        FQ_dict = 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.qc_hybrid_bwa_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        SAM = 'qc_hybrid/{name}/read_mapping/illumina/bwa_readmap_sam.io',
        INFORMS = 'qc_hybrid/{name}/read_mapping/illumina/commands.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    threads: 8
    run:
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput
        from camel.app.tools.bwa.bwamap import BWAMap
        bwa_map = BWAMap()
        fq_in = FastqInput.from_fq_dict(Path(input.FQ_dict), 'illumina')
        bwa_map.add_input_files({'FASTQ_PE': fq_in.pe})
        bwa_map.update_parameters(threads=threads)
        snakemakeutils.add_pickle_input(bwa_map, 'INDEX_GENOME_PREFIX', Path(input.INDEX_GENOME_PREFIX_BWA))
        step = Step(str(rule), bwa_map, dir_=Path(str(params.dir_)))
        step.run()
        bwa_map.informs['_tag'] = wildcards.name
        snakemakeutils.dump_tool_outputs(bwa_map, output)

rule qc_hybrid_read_mapping_ont:
    """
    Maps the long reads against the assembly.
    """
    input:
        FQ = 'fq_dict.io',
        INDEX_GENOME_PREFIX_BWA = rules.qc_hybrid_bwa_index.output.INDEX_GENOME_PREFIX,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        SAM = 'qc_hybrid/{name}/read_mapping/ont/ont_mapping.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    threads: 8
    run:
        from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping
        minimap2 = Minimap2Mapping()
        snakemakeutils.add_pickle_input(minimap2, 'FASTA', Path(input.FASTA))
        minimap2.add_input_files(snakepipelineutils.extract_fq_input(Path(input.FQ), key_se='FASTQ', read_type='SE'))
        step = Step(str(rule), minimap2, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(minimap2, output)

rule qc_hybrid_sam_to_bam_illumina:
    """
    Converts SAM to BAM from the short reads mapping.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_illumina.output.SAM
    output:
        BAM =  'qc_hybrid/{name}/read_mapping/illumina/bwa_readmap_bam.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView()
        snakemakeutils.add_pickle_inputs(samtools_view, input)
        step = Step(str(rule), samtools_view, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_view, output)

rule qc_hybrid_bam_sorting_illumina:
    """
    Sorts the BAM alignment from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_sam_to_bam_illumina.output.BAM
    output:
        BAM = 'qc_hybrid/{name}/read_mapping/illumina/bwa_readmap_sorted.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort()
        snakemakeutils.add_pickle_inputs(samtools_sort, input)
        step = Step(str(rule), samtools_sort, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_sort, output)

rule qc_hybrid_bam_indexing_illumina:
    """
    Indexes the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM
    output:
        BAM = 'qc_hybrid/{name}/read_mapping/illumina/samtools-index.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex()
        snakemakeutils.add_pickle_inputs(samtools_index, input)
        step = Step(str(rule), samtools_index, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_index, output)

rule qc_hybrid_mapping_stats_illumina:
    """
    Retrieves mapping stats from the short read mappings.
    """
    input:
        BAM = rules.qc_hybrid_bam_indexing_illumina.output.BAM
    output:
        INFORMS = 'qc_hybrid/{name}/read_mapping/illumina/flagstat.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        samtools_flagstat = SamtoolsFlagstat()
        snakemakeutils.add_pickle_inputs(samtools_flagstat, input)
        step = Step(str(rule), samtools_flagstat, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_flagstat, output)

rule qc_hybrid_samtools_depth_illumina:
    """
    Runs samtools depth on the BAM file from the short reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_illumina.output.BAM
    output:
        TSV = 'qc_hybrid/{name}/read_mapping/illumina/tsv.io',
        INFORMS = 'qc_hybrid/{name}/read_mapping/illumina/samtools-depth.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/illumina'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth()
        snakemakeutils.add_pickle_inputs(samtools_depth, input)
        step = Step(str(rule), samtools_depth, dir_=Path(str(params.dir_)))
        step.run()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        snakemakeutils.dump_tool_outputs(samtools_depth, output)

rule qc_hybrid_sam_to_bam_ont:
    """
    Converts SAM to BAM from the long reads mapping.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_ont.output.SAM
    output:
        BAM = 'qc_hybrid/{name}/read_mapping/ont/minimap2_bam.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    run:
        from camel.app.tools.samtools.samtoolsview import SamtoolsView
        samtools_view = SamtoolsView()
        snakemakeutils.add_pickle_inputs(samtools_view, input)
        step = Step(str(rule), samtools_view, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_view, output)

rule qc_hybrid_bam_sorting_ont:
    """
    Sorts the BAM alignment from the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_sam_to_bam_ont.output.BAM
    output:
        BAM = 'qc_hybrid/{name}/read_mapping/ont/minimap2_sorted.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    run:
        from camel.app.tools.samtools.samtoolssort import SamtoolsSort
        samtools_sort = SamtoolsSort()
        snakemakeutils.add_pickle_inputs(samtools_sort, input)
        step = Step(str(rule), samtools_sort, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_sort, output)

rule qc_hybrid_bam_indexing_ont:
    """
    Indexes the BAM file from the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM
    output:
        BAM = 'qc_hybrid/{name}/read_mapping/ont/bwa-index-longreads.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    run:
        from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
        samtools_index = SamtoolsIndex()
        snakemakeutils.add_pickle_inputs(samtools_index, input)
        step = Step(str(rule), samtools_index, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_index, output)

rule qc_hybrid_mapping_stats_ont:
    """
    Retrieves mapping stats from the long read mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM
    output:
        INFORMS = 'qc_hybrid/{name}/read_mapping/ont/flagstat-longreads.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    run:
        from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
        samtools_flagstat = SamtoolsFlagstat()
        snakemakeutils.add_pickle_inputs(samtools_flagstat, input)
        step = Step(str(rule), samtools_flagstat, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(samtools_flagstat, output)

rule qc_hybrid_samtools_depth_ont:
    """
    Runs samtools depth on the BAM file of the long reads mapping.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM
    output:
        TSV = 'qc_hybrid/{name}/read_mapping/ont/tsv-long.io',
        INFORMS = 'qc_hybrid/{name}/read_mapping/ont/samtools-depth-long.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/read_mapping/ont'
    run:
        from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
        samtools_depth = SamtoolsDepth()
        snakemakeutils.add_pickle_inputs(samtools_depth, input)
        step = Step(str(rule), samtools_depth, dir_=Path(str(params.dir_)))
        step.run()
        samtools_depth.informs['_tag'] = 'Coverage calculation'
        snakemakeutils.dump_tool_outputs(samtools_depth, output)

rule qc_hybrid_freebayes:
    """
    Checks for small variants in the assembly.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_illumina.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_illumina.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF =  'qc_hybrid/{name}/freebayes/vcf.io',
        INFORMS = 'qc_hybrid/{name}/freebayes/commands.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/freebayes',
        freebayes_options = config.get('freebayes',{})
    run:
        from camel.app.tools.freebayes.freebayes import Freebayes
        freebayes = Freebayes()
        snakemakeutils.add_pickle_inputs(freebayes, input)
        freebayes.update_parameters(**params.freebayes_options)
        step = Step(str(rule), freebayes, dir_=Path(str(params.dir_)))
        step.run()
        freebayes.informs['_tag'] = wildcards.name
        snakemakeutils.dump_tool_outputs(freebayes, output)

rule qc_hybrid_sniffles:
    """
    Checks for structural variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF = 'qc_hybrid/{name}/sniffles/vcf.io',
        INFORMS = 'qc_hybrid/{name}/sniffles/commands.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/sniffles'
    threads: 4
    run:
        from camel.app.tools.sniffles.sniffles import Sniffles
        sniffles = Sniffles()
        snakemakeutils.add_pickle_inputs(sniffles, input)
        sniffles.update_parameters(threads=threads)
        step = Step(str(rule), sniffles, dir_=Path(str(params.dir_)))
        step.run()
        sniffles.informs['_tag'] = wildcards.name
        snakemakeutils.dump_tool_outputs(sniffles, output)

rule qc_hybrid_clair3:
    """
    Checks for small variants in the assembly using long reads.
    """
    input:
        BAM = rules.qc_hybrid_bam_sorting_ont.output.BAM,
        BAM_INDEX = rules.qc_hybrid_bam_indexing_ont.output.BAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        VCF = 'qc_hybrid/{name}/clair3_output/gzipped_vcf.io',
        INFORMS = 'qc_hybrid/{name}/clair3_output/commands.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}',
        clair3_options = config.get('clair3', {})
    threads: 8
    run:
        from camel.app.tools.clair3.clair3 import Clair3
        clair3 = Clair3()
        snakemakeutils.add_pickle_inputs(clair3, input)
        clair3.update_parameters(
            **params.clair3_options, chunk_size=100_000, platform='ont', no_phasing=True, include_ctgs=True,
            threads=threads)
        step = Step(str(rule), clair3, dir_=Path(str(params.dir_)))
        step.run()
        clair3.informs['_tag'] = wildcards.name
        snakemakeutils.dump_tool_outputs(clair3, output)

rule qc_hybrid_unzip_clair3_vcf:
    """
    Unzips the output from clair3 to be compatible with the qc_add_vcf_info_to_informs rule.
    """
    input:
        VCF_GZ = rules.qc_hybrid_clair3.output.VCF
    output:
        VCF = 'qc_hybrid/{name}/clair3_output/vcf.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/clair3_output'
    run:
        from camel.app.tools.bcftools.bcftoolsview import BcftoolsView
        bcftools_view = BcftoolsView()
        snakemakeutils.add_pickle_inputs(bcftools_view, input)
        bcftools_view.update_parameters(compress_output=False, output_filename='variants.vcf')
        step = Step(str(rule), bcftools_view, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(bcftools_view, output)

rule qc_hybrid_add_vcf_info_to_informs:
    """
    Parses the VCF file and adds the statistics to the informs.
    """
    input:
        INFORMS = 'qc_hybrid/{name}/{method}/commands.io',
        VCF = 'qc_hybrid/{name}/{method}/vcf.io'
    output:
        INFORMS = 'qc_hybrid/{name}/{method}/informs.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/{wildcards.method}'
    run:
        from camel.app.core.utils import vcfutils
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        if wildcards.method == 'sniffles':
            informs['nb_of_variants'] = sum(informs['variants'].values())
            informs['nb_of_snps'] = 0
            informs['nb_of_indels'] = informs['variants']['DEL'] + informs['variants']['INS']
            informs['nb_of_svs'] = informs['variants']['BND'] + informs['variants']['DUP'] + informs['variants']['INV']
        else:
            variants = vcfutils.retrieve_variants(snakemakeutils.load_object(Path(input.VCF))[0].path)
            informs['nb_of_variants'] = len(variants)
            informs['nb_of_snps'] = sum(v.is_snp for v in variants)
            informs['nb_of_indels'] = sum(v.is_indel for v in variants)
            informs['nb_of_svs'] = 0
        snakemakeutils.dump_object(informs, Path(output.INFORMS))

rule qc_hybrid_ale:
    """
    Generates ALE QC report for the final assembly.
    """
    input:
        SAM = rules.qc_hybrid_read_mapping_illumina.output.SAM,
        FASTA = rules.qc_hybrid_samtools_index.output.FASTA
    output:
        ALE = 'qc_hybrid/{name}/ale_illumina/ALE.io',
        INFORMS = 'qc_hybrid/{name}/ale_illumina/informs-report.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/ale_illumina'
    run:
        from camel.app.tools.ale.ale import ALE
        ale_report = ALE()
        snakemakeutils.add_pickle_inputs(ale_report, input)
        step = Step(str(rule), ale_report, dir_=Path(str(params.dir_)))
        step.run()
        ale_report.informs['_tag'] = f'{wildcards.name}'
        snakemakeutils.dump_tool_outputs(ale_report, output)

rule qc_hybrid_ale2wiggle:
    """
    Generates wiggle files from ALE for IGV visualization.
    """
    input:
        ALE = rules.qc_hybrid_ale.output.ALE
    output:
        TSV = 'qc_hybrid/{name}/ale_illumina/wiggle.io',
        INFORMS = 'qc_hybrid/{name}/ale_illumina/commands-wiggle.io'
    params:
        dir_ = lambda wildcards: f'qc_hybrid/{wildcards.name}/ale_illumina'
    run:
        from camel.app.tools.ale.ale2wiggle import ALE2Wiggle
        ale2wiggle_report = ALE2Wiggle()
        snakemakeutils.add_pickle_inputs(ale2wiggle_report, input)
        step = Step(str(rule), ale2wiggle_report, dir_=Path(str(params.dir_)))
        step.run()
        ale2wiggle_report.informs['_tag'] = f'{wildcards.name}'
        snakemakeutils.dump_tool_outputs(ale2wiggle_report, output)
