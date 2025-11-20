import json
import logging
import shutil
from pathlib import Path
from typing import Union

import pandas as pd

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping


def _get_reference_genome(wildcards) -> str:
    """
    Returns the reference genome used for the current iteration of the iterative mapping.
    :param wildcards: Snakemake wildcards
    :return: Path to reference genome IO file
    """
    iteration_nb = int(wildcards.nb_iter)
    if iteration_nb == 1:
        # The first iteration uses the provided reference genome
        return iterativemapping.INPUT_FASTA_REF
    else:
        # The following iterations map to the consensus sequence generated in the previous iteration
        iteration_nb_prev = f'{iteration_nb-1:02d}'
        logging.info(f'Using FASTA file from previous iteration as reference genome ({iteration_nb_prev})')
        return f'iterative_mapping/iter_{iteration_nb_prev}/phase_2-variants/fasta.io'

rule iterative_mapping_phase_1_map_reads:
    """
    Maps the reads against the input sequence.
    """
    input:
        FASTA = lambda wildcards: _get_reference_genome(wildcards),
        IO_FASTQ = 'preprocess/downsample/fq_dict.io'
    output:
        BAM = 'iterative_mapping/iter_{nb_iter}/phase_1-mapping/bam.io',
        JSON = 'iterative_mapping/iter_{nb_iter}/phase_1-mapping/stats.json',
        INFORMS = 'iterative_mapping/iter_{nb_iter}/phase_1-mapping/informs.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_1-mapping',
        input_type = config['input']['type'],
        gap_depth_cutoff = config['low_depth'].get('gap_depth_cutoff'),
        gap_len_cutoff = config['low_depth'].get('gap_len_cutoff')
    threads: 16
    run:
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow

        # Run the workflow
        workflow = ReadMappingWorkflow(Path(str(params.dir_)).absolute())
        out = workflow.run(
            fastq_in=FastqInput.from_fq_dict(Path(input.IO_FASTQ), params.input_type),
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            gap_len_cutoff=params.gap_len_cutoff,
            gap_depth_cutoff=params.gap_depth_cutoff,
            threads=threads)

        # Save output
        snakemakeutils.dump_object([ToolIOFile(out.path_bam)], Path(output.BAM))
        with open(output.JSON, 'w') as handle:
            json.dump(out.stats, handle, indent=2)
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule iterative_mapping_phase_1_call_variants_bcftools:
    """
    Calls variants using bcftools.
    """
    input:
        FASTA = lambda wildcards: _get_reference_genome(wildcards),
        BAM = lambda wildcards: rules.iterative_mapping_phase_1_map_reads.output.BAM
    output:
        VCF = 'iterative_mapping/iter_{nb_iter}/phase_1-variants/vcf.io',
        INFORMS = 'iterative_mapping/iter_{nb_iter}/phase_1-variants/informs-call.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_1-variants',
        input_type = config['input']['type']
    threads: 1
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.callvariants import CallVariants
        workflow = CallVariants(Path(str(params.dir_)).absolute())
        out = workflow.run(
            bam_in=snakemakeutils.load_object(Path(str(input.BAM)))[0].path,
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            input_type=params.input_type,
            caller='bcftools',
            params={},
            threads=threads
        )
        snakemakeutils.dump_object([ToolIOFile(out.path_vcf)], Path(output.VCF))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule iterative_mapping_phase_1_filter_variants_bcftools:
    """
    Filters the variants called by bcftools.
    """
    input:
        VCF = rules.iterative_mapping_phase_1_call_variants_bcftools.output.VCF
    output:
        VCF = 'iterative_mapping/iter_{nb_iter}/phase_1-variants/vcf-filt.io',
        INFORMS = 'iterative_mapping/iter_{nb_iter}/phase_1-variants/informs-filt.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_1-variants',
        min_af = config['iterative_mapping'].get('variant_filters', {}).get('min_af', 0.5),
        min_dp = config['iterative_mapping'].get('variant_filters', {}).get('min_dp', 5),
        min_qual = config['iterative_mapping'].get('variant_filters', {}).get('min_qual', 25),
        min_mq = config['iterative_mapping'].get('variant_filters', {}).get('min_mq', 30)
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.filtervariants import FilterVariants
        workflow = FilterVariants(Path(str(params.dir_)).absolute())
        out = workflow.run(
            vcf_in=snakemakeutils.load_object(Path(str(input.VCF)))[0].path,
            calling_method='bcftools',
            filters={
                'min_af': params.min_af,
                'min_dp': params.min_dp,
                'min_qual': params.min_qual,
                'min_mq': params.min_mq})
        snakemakeutils.dump_object([ToolIOFile(out.path_vcf)], Path(output.VCF))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule iterative_mapping_phase_1_apply_variants:
    """
    Applies the phase 1 variants to obtain an updated consensus sequence.
    """
    input:
        FASTA = lambda wildcards: _get_reference_genome(wildcards),
        VCF = rules.iterative_mapping_phase_1_filter_variants_bcftools.output.VCF
    output:
        FASTA = 'iterative_mapping/iter_{nb_iter}/phase_1-variants/fasta.io'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_1-variants',
        name = lambda wildcards: config['input']['sample_name'],
        iter_nb = lambda wildcards: wildcards.nb_iter
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.applyvariants import ApplyVariants

        apply_variants = ApplyVariants(Path(str(params.dir_)).absolute())
        out = apply_variants.run(
            fasta_in=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            vcf_in=snakemakeutils.load_object(Path(str(input.VCF)))[0].path,
            name=str(params.name),
            description=f'iter_{params.iter_nb}-phase_1'
        )
        snakemakeutils.dump_object([ToolIOFile(out.path_fasta)], Path(output.FASTA))

rule iterative_mapping_phase_2_map_reads:
    """
    Maps the reads against the consensus sequence generated in phase 1.
    """
    input:
        FASTA = rules.iterative_mapping_phase_1_apply_variants.output.FASTA,
        IO_FASTQ = 'preprocess/downsample/fq_dict.io'
    output:
        BAM = 'iterative_mapping/iter_{nb_iter}/phase_2-mapping/bam.io',
        BED = 'iterative_mapping/iter_{nb_iter}/phase_2-mapping/bed.io',
        JSON = 'iterative_mapping/iter_{nb_iter}/phase_2-mapping/informs.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_2-mapping',
        input_type = config['input']['type'],
        gap_depth_cutoff = config['low_depth'].get('gap_depth_cutoff'),
        gap_len_cutoff = config['low_depth'].get('gap_len_cutoff')
    threads: 16
    run:
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow

        # Run the workflow
        workflow = ReadMappingWorkflow(Path(str(params.dir_)).absolute())
        out = workflow.run(
            fastq_in=FastqInput.from_fq_dict(Path(input.IO_FASTQ), params.input_type),
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            gap_len_cutoff=params.gap_len_cutoff,
            gap_depth_cutoff=params.gap_depth_cutoff,
            threads=threads)

        # Save output
        snakemakeutils.dump_object([ToolIOFile(out.path_bam)], Path(output.BAM))
        snakemakeutils.dump_object([ToolIOFile(out.path_bed_low_cov)], Path(output.BED))
        with open(output.JSON, 'w') as handle:
            json.dump(out.stats, handle, indent=2)

rule iterative_mapping_phase_2_call_variants_clair3:
    """
    Calls variants using Clair3.
    """
    input:
        FASTA = rules.iterative_mapping_phase_1_apply_variants.output.FASTA,
        BAM = lambda wildcards: rules.iterative_mapping_phase_2_map_reads.output.BAM
    output:
        VCF = 'iterative_mapping/iter_{nb_iter}/phase_2-variants/vcf.io',
        INFORMS= 'iterative_mapping/iter_{nb_iter}/phase_2-variants/informs-call.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_2-variants',
        input_type = config['input']['type'],
        model_path = config['iterative_mapping'].get('clair3', {}).get('model')
    threads: 4
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.callvariants import CallVariants

        # Run workflow
        workflow = CallVariants(Path(str(params.dir_)).absolute())
        out = workflow.run(
            bam_in=snakemakeutils.load_object(Path(str(input.BAM)))[0].path,
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            input_type=params.input_type,
            caller='clair3',
            params={'model': Path(params.model_path)},
            threads=threads
        )
        snakemakeutils.dump_object([ToolIOFile(out.path_vcf)],Path(output.VCF))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule iterative_mapping_phase_2_filter_variants_clair3:
    """
    Filters the variants called by Clair3.
    """
    input:
        VCF = rules.iterative_mapping_phase_2_call_variants_clair3.output.VCF
    output:
        VCF = 'iterative_mapping/iter_{nb_iter}/phase_2-variants/vcf-filt.io',
        INFORMS = 'iterative_mapping/iter_{nb_iter}/phase_2-variants/informs-filt.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_2-variants',
        min_af = config['iterative_mapping'].get('variant_filters', {}).get('min_af', 0.5),
        min_dp = config['iterative_mapping'].get('variant_filters', {}).get('min_dp', 5),
        min_qual = config['iterative_mapping'].get('variant_filters', {}).get('min_qual', 25)
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.filtervariants import FilterVariants
        workflow = FilterVariants(Path(str(params.dir_)).absolute())
        out = workflow.run(
            vcf_in=snakemakeutils.load_object(Path(str(input.VCF)))[0].path,
            calling_method='clair3',
            filters={'min_af': params.min_af, 'min_dp': params.min_dp, 'min_qual': params.min_qual})
        snakemakeutils.dump_object([ToolIOFile(out.path_vcf)], Path(output.VCF))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule iterative_mapping_phase_2_apply_variants:
    """
    Applies the phase 2 variants to obtain the final consensus sequence for this iteration.
    """
    input:
        FASTA = rules.iterative_mapping_phase_1_apply_variants.output.FASTA,
        VCF = rules.iterative_mapping_phase_2_filter_variants_clair3.output.VCF
    output:
        FASTA = 'iterative_mapping/iter_{nb_iter}/phase_2-variants/fasta.io',
        INFORMS = 'iterative_mapping/iter_{nb_iter}/phase_2-variants/informs-apply.json'
    params:
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}/phase_2-variants',
        name = lambda wildcards: config['input']['sample_name'],
        iter_nb = lambda wildcards: wildcards.nb_iter
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.applyvariants import ApplyVariants
        apply_variants = ApplyVariants(Path(str(params.dir_)).absolute())
        out = apply_variants.run(
            fasta_in=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            vcf_in=snakemakeutils.load_object(Path(str(input.VCF)))[0].path,
            name=str(params.name),
            description=f'iter_{params.iter_nb}-phase_2'
        )
        snakemakeutils.dump_object([ToolIOFile(out.path_fasta)], Path(output.FASTA))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

checkpoint iterative_mapping_collect_stats:
    """
    Collects the stats for a single iteration in the iterative mapping. 
    """
    input:
        VCF_filt_p1 = rules.iterative_mapping_phase_1_filter_variants_bcftools.output.VCF,
        VCF_filt_p2 = rules.iterative_mapping_phase_2_filter_variants_clair3.output.VCF,
        FASTA = rules.iterative_mapping_phase_2_apply_variants.output.FASTA,
        JSON = rules.iterative_mapping_phase_2_map_reads.output.JSON
    output:
        JSON = 'iterative_mapping/iter_{nb_iter}/stats-iter_{nb_iter}.json'
    params:
        nb_iter = lambda wildcards: wildcards.nb_iter,
        dir_ = lambda wildcards: f'iterative_mapping/iter_{wildcards.nb_iter}'
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.viral_consensus.combineiterativemappingstats import CollectIterativeMappingStats
        combiner = CollectIterativeMappingStats()
        combiner.add_input_files({
            'FASTA': snakemakeutils.load_object(Path(input.FASTA)),
            'JSON_depth': [ToolIOFile(Path(input.JSON))],
            'VCF_p1': snakemakeutils.load_object(Path(input.VCF_filt_p1)),
            'VCF_p2': snakemakeutils.load_object(Path(input.VCF_filt_p2))
        })
        step = Step(rule_name=str(rule), tool=combiner, dir_=Path(str(params.dir_)))
        combiner.update_parameters(nb_iter=str(params.nb_iter), output_filename=str(Path(output.JSON).name))
        step.run()
        shutil.copyfile(combiner.tool_outputs['JSON'][0].path, output.JSON)

def _check_if_converged(wildcards) -> Union[Path, list[Path]]:
    """
    Checks if the consensus sequence converged (no changes in the last two iterations).
    :return: Path / list of paths to generate
    """
    # Get the directory of the latest iteration
    dirs_iter = sorted(list(Path('iterative_mapping').glob('iter_*')))

    # Return the directory for the first iteration if there are no directories yet
    if len(dirs_iter) == 0:
        return Path('iterative_mapping/iter_01/stats-iter_01.json')

    # Return the directory for the second iteration after the first iteration
    if len(dirs_iter) == 1:
        return Path('iterative_mapping/iter_02/stats-iter_02.json')

    # Get the sequence hash for each iteration
    hash_by_iter = {}
    for path_json in sorted([dir_ / f'stats-{dir_.name}.json' for dir_ in dirs_iter]):
        with path_json.open() as h:
            data_stats = json.load(h)
        hash_by_iter[path_json.parent.name] = data_stats['all_segments']['sequence_md5']
    logging.info(f'Hashes: {hash_by_iter}')

    # Check if the number of variants == 0 for the last two iterations
    hashes_last_two_iter = [hash_by_iter[dirs_iter[i].name] for i in (-1, -2)]
    has_converged = len(set(hashes_last_two_iter)) == 1
    if has_converged or len(dirs_iter) >= config['iterative_mapping'].get('max_iter', 6):
        return [dir_ / f'stats-{dir_.name}.json' for dir_ in dirs_iter]

    # Otherwise, return the next iteration
    iter_nb = int(dirs_iter[-1].name.split('_')[-1])
    return Path('iterative_mapping', f'iter_{iter_nb+1:02d}', f'stats-iter_{iter_nb+1:02d}.json')

rule iterative_mapping_combine_stats_all_iterations:
    """
    Combines the statistics for each iteration of the mapping approach.
    """
    input:
        JSON = _check_if_converged
    output:
        TSV = 'iterative_mapping/stats/stats_all_iterations.tsv',
        TSV_seg = 'iterative_mapping/stats/stats_all_iterations-by_segment.tsv'
    run:
        # Combined across all segments
        records_out = []
        records_out_by_segment = []
        # noinspection PyTypeChecker
        for path_json in [Path(x) for x in input.JSON]:
            nb_iter = int(path_json.parent.name.replace('iter_', ''))

            # Global
            with path_json.open() as handle:
                data = json.load(handle)
            records_out.append(data['all_segments'])

            # By segment
            for seq_id, data_segment in data['by_segment'].items():
                records_out_by_segment.append({
                    'seq_id': seq_id,
                    'segment': seq_id.split('-')[-1],
                    'iter': nb_iter,
                    **data_segment
                })
        pd.DataFrame(records_out).to_csv(output.TSV, sep='\t', index=False)
        pd.DataFrame(records_out_by_segment).to_csv(output.TSV_seg, sep='\t', index=False)

rule iterative_mapping_select_output:
    """
    Selects the output of the iterative mapping workflow.
    """
    input:
        TSV = rules.iterative_mapping_combine_stats_all_iterations.output.TSV
    output:
        BAM = 'iterative_mapping/output/bam.io',
        BED = 'iterative_mapping/output/bed.io',
        VCF_p1 = 'iterative_mapping/output/vcf_p1.io',
        VCF_p2 = 'iterative_mapping/output/vcf_p2.io',
        FASTA = 'iterative_mapping/output/fasta.io',
        FASTA_ref = 'iterative_mapping/output/fasta-ref.io'
    run:
        data_mapping = pd.read_table(input.TSV, dtype=str)

        # Consensus FASTA file
        path_fasta = str(rules.iterative_mapping_phase_2_apply_variants.output.FASTA).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_fasta, output.FASTA)

        # Reference FASTA file
        path_fasta = str(rules.iterative_mapping_phase_1_apply_variants.output.FASTA).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_fasta, output.FASTA_ref)

        # BAM file
        path_bam = str(rules.iterative_mapping_phase_2_map_reads.output.BAM).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_bam,output.BAM)

        # Low depth BED file
        path_bed = str(rules.iterative_mapping_phase_2_map_reads.output.BED).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_bed, output.BED)

        # VCF file (P1)
        path_vcf = str(rules.iterative_mapping_phase_1_filter_variants_bcftools.output.VCF).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_vcf, output.VCF_p1)

        # VCF file (P2)
        path_vcf = str(rules.iterative_mapping_phase_2_filter_variants_clair3.output.VCF).format(
            nb_iter=data_mapping['iter'].iloc[-1])
        shutil.copyfile(path_vcf, output.VCF_p2)

rule iterative_mapping_trim_fasta_edges:
    """
    Trims the edges of the FASTA file where sequencing depth is too low.
    """
    input:
        FASTA = rules.iterative_mapping_select_output.output.FASTA,
        BED = rules.iterative_mapping_select_output.output.BED
    output:
        FASTA = 'iterative_mapping/trim_edges/fasta-trim.io',
        INFORMS = 'iterative_mapping/trim_edges/informs.io'
    params:
        dir_ = 'iterative_mapping/output'
    run:
        from camel.app.tools.pipelines.viral_consensus.trimfastaedges import TrimFastaEdges
        trim_edges = TrimFastaEdges()
        snakemakeutils.add_pickle_inputs(trim_edges, input)
        step = Step(rule_name=str(rule), tool=trim_edges, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(trim_edges, output)

rule iterative_mapping_report:
    """
    Generates the HTML report for the iterative mapping workflow.
    """
    input:
        FASTA = rules.iterative_mapping_trim_fasta_edges.output.FASTA,
        FASTA_ref = rules.iterative_mapping_select_output.output.FASTA_ref,
        BAM = rules.iterative_mapping_select_output.output.BAM,
        TSV = rules.iterative_mapping_combine_stats_all_iterations.output.TSV,
        TSV_seg = rules.iterative_mapping_combine_stats_all_iterations.output.TSV_seg,
        VCF_p1 = rules.iterative_mapping_select_output.output.VCF_p1,
        VCF_p2 = rules.iterative_mapping_select_output.output.VCF_p2,
        INFORMS_trim_fasta = rules.iterative_mapping_trim_fasta_edges.output.INFORMS
    output:
        VAL_HTML = 'iterative_mapping/report/html.iob' # iterativemapping.OUTPUT_REPORT
    params:
        name = config['input']['sample_name'],
        dir_ = 'iterative_mapping/report',
        max_iter = config['iterative_mapping'].get('max_iter', 6),
        gap_depth_cutoff = config['low_depth'].get('gap_depth_cutoff', 50),
        gap_len_cutoff = config['low_depth'].get('gap_len_cutoff', 10)
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.viral_consensus.reporteriterativemapping import ReporterIterativeMapping
        reporter = ReporterIterativeMapping()
        reporter.add_input_files({
            'TSV': [ToolIOFile(Path(input.TSV))],
            'TSV_seg': [ToolIOFile(Path(input.TSV_seg))]
        })
        snakemakeutils.add_pickle_inputs(reporter, input,
            ['FASTA', 'FASTA_ref', 'BAM', 'VCF_p1', 'VCF_p2', 'INFORMS_trim_fasta'])
        reporter.update_parameters(
            name=str(params.name), max_iter=params.max_iter, gap_depth_cutoff=params.gap_depth_cutoff,
            gap_len_cutoff=params.gap_len_cutoff)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule iterative_mapping_summary:
    """
    Creates the summary output for the iterative mapping workflow.
    """
    input:
        TSV = 'iterative_mapping/stats/stats_all_iterations.tsv',
        INFORMS_low_depth = rules.iterative_mapping_trim_fasta_edges.output.INFORMS
    output:
        FILE = 'iterative_mapping/summary/iterative_mapping.{ext}' # iterativemapping.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        data_in = pd.read_table(input.TSV)

        # Low depth
        informs_depth = snakemakeutils.load_object(Path(input.INFORMS_low_depth))

        # Iterative mapping
        rows_out = [
            ('iterative_mapping_nb_iterations', len(data_in)),
            ('iterative_mapping_sequence_length', int(data_in.iloc[-1]['length'])),
            ('iterative_mapping_median_depth', data_in.iloc[-1]['depth_median']),
            ('iterative_mapping_median_depth_iqr', data_in.iloc[-1]['depth_iqr']),
            ('iterative_mapping_covered_rate', data_in.iloc[-1]['covered_rate']),
            ('low_depth_nb_clipped', int(informs_depth['nb_clipped'])),
            ('low_depth_nb_masked', int(informs_depth['nb_masked']))
        ]

        # Save to summary
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'preprocess')

rule iterative_mapping_combine_informs:
    """
    Combines the informs for the iterative mapping workflow.
    """
    input:
        INFORMS_mapping = str(rules.iterative_mapping_phase_1_map_reads.output.INFORMS).format(nb_iter='01'),
        INFORMS_calling_p1 = str(rules.iterative_mapping_phase_1_call_variants_bcftools.output.INFORMS).format(nb_iter='01'),
        INFORMS_calling_p2 = str(rules.iterative_mapping_phase_2_call_variants_clair3.output.INFORMS).format(nb_iter='01'),
        INFORMS_filtering_p1 = str(rules.iterative_mapping_phase_1_filter_variants_bcftools.output.INFORMS).format(nb_iter='01'),
        INFORMS_filtering_p2 = str(rules.iterative_mapping_phase_2_filter_variants_clair3.output.INFORMS).format(nb_iter='01'),
        INFORMS_apply = str(rules.iterative_mapping_phase_2_apply_variants.output.INFORMS).format(nb_iter='01')
    output:
        INFORMS = 'iterative_mapping/report/informs.io' # iterativemapping.OUTPUT_INFORMS
    run:
        data_out = []
        # noinspection PyUnresolvedReferences
        for key, path_json in input.items():

            # Add suffix
            if str(key).endswith('_p1'):
                suffix = ' - Phase 1 (bcftools)'
            elif str(key).endswith('_p2'):
                suffix = ' - Phase 2 (Clair3)'
            else:
                suffix = None

            # Save output
            with Path(path_json).open() as handle:
                informs_all = json.load(handle)
                for inform in informs_all:
                    inform['_tag'] = 'Iterative mapping' + (suffix if suffix is not None else '')
                    data_out.append(inform)
        snakemakeutils.dump_object(data_out, Path(output.INFORMS))
