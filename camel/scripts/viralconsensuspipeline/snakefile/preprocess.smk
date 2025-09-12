import json
from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.scripts.viralconsensuspipeline.snakefile import iterativemapping


#############
# AmpliGone #
#############
rule preprocess_ampligone_to_bed:
    """
    Creates a BED file with the location of the primers in the reference genome.
    """
    input:
        FASTA_ref = iterativemapping.INPUT_FASTA_REF,
        FASTA_primers = config.get('preprocess', {}).get('ampligone', {}).get('fasta', [])
    output:
        BED = 'preprocess/ampligone/tool/bed.io',
        INFORMS = 'preprocess/ampligone/tool/informs.io'
    params:
        dir_ = 'preprocess/ampligone/tool',
        primer_mismatch_rate = '0.01' if config['input_type'] == 'illumina' else '0.1'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.ampligone.ampligonefasta2bed import AmpliGoneFasta2Bed
        ampligone = AmpliGoneFasta2Bed()
        snakemakeutils.add_pickle_input(ampligone, 'FASTA_ref', Path(input.FASTA_ref))
        if not Path(input.FASTA_primers).exists():
            raise FileNotFoundError(f'Cannot find FASTA file with primers: {input.FASTA_primers}')
        ampligone.add_input_files({'FASTA_primers': [ToolIOFile(Path(input.FASTA_primers))]})
        step = Step(rule_name=str(rule), tool=ampligone, dir_=Path(str(params.dir_)))
        ampligone.update_parameters(primer_mismatch_rate=params.primer_mismatch_rate)
        step.run()
        snakemakeutils.dump_tool_outputs(ampligone, output)

rule preprocess_remove_ampligone_report:
    """
    Creates a report for the primer removal step.
    """
    input:
        BED = rules.preprocess_ampligone_to_bed.output.BED,
        INFORMS_ampligone = rules.preprocess_ampligone_to_bed.output.INFORMS
    output:
        HTML = 'preprocess/ampligone/report/html.iob' # OUTPUT_AMPLIGONE_REPORT
    params:
        dir_ = 'preprocess/ampligone/report'
    threads: 4
    run:
        from camel.app.tools.ampligone.ampligonefasta2bedreporter import AmpliGoneFasta2BedReporter
        reporter = AmpliGoneFasta2BedReporter()
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule preprocess_remove_ampligone_report_empty:
    """
    Creates an empty output report for the AmpliGone step.
    """
    output:
        HTML = 'preprocess/ampligone/report/html-empty.iob' # preprocess.OUTPUT_AMPLIGONE_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Primer localization', Path(output.HTML))

############################
# Per-segment downsampling #
############################
rule preprocess_map_reads_pre:
    """
    Maps the reads against the reference genome.
    """
    input:
        FASTA = iterativemapping.INPUT_FASTA_REF,
        IO_FASTQ = 'fq_dict.io'
    output:
        BAM = 'preprocess/mapping_pre/bam.io',
        JSON = 'preprocess/mapping_pre/stats.json',
        INFORMS = 'preprocess/mapping_pre/informs.json'
    params:
        dir_ = lambda wildcards: 'preprocess/mapping_pre',
        input_type = config['input_type'],
        gap_depth_cutoff = config['low_depth'].get('gap_depth_cutoff', 5),
        gap_len_cutoff = config['low_depth'].get('gap_len_cutoff', 10)
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow

        # Run the workflow
        workflow = ReadMappingWorkflow(Path(str(params.dir_)).absolute())
        out = workflow.run(
            fastq_in=FastqInput.from_fq_dict(Path(input.IO_FASTQ), params.input_type),
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            threads=threads,
            gap_depth_cutoff=int(params.gap_depth_cutoff),
            gap_len_cutoff=int(params.gap_len_cutoff),
        )

        # Save output
        snakemakeutils.dump_object([ToolIOFile(out.path_bam)], Path(output.BAM))
        with open(output.JSON, 'w') as handle:
            json.dump(out.stats, handle, indent=2)
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule preprocess_downsample_by_segment:
    """
    Down-samples the FASTQ files by segment. 
    """
    input:
        JSON = rules.preprocess_map_reads_pre.output.JSON,
        BAM = rules.preprocess_map_reads_pre.output.BAM,
        BED = rules.preprocess_ampligone_to_bed.output.BED if 'ampligone' in config['analyses'] else []
    output:
        FASTQ = 'preprocess/downsample/fq_dict.io',
        INFORMS = 'preprocess/downsample/informs.json'
    params:
        dir_ = 'preprocess/downsample',
        input_type = config['input_type'],
        max_depth = config['downsampling'].get('coverage_max_by_segment', 250)
    threads: 8
    run:
        from camel.scripts.viralconsensuspipeline.workflows.segmentdownsampling import SegmentDownsamplingWorkflow
        workflow = SegmentDownsamplingWorkflow(Path(params.dir_).absolute())
        path_bam = snakemakeutils.load_object(Path(input.BAM))[0].path
        path_bed = snakemakeutils.load_object(Path(input.BED))[0].path if len(input.BED) > 0 else None
        out = workflow.run(path_bam, params.input_type, Path(input.JSON), params.max_depth, path_bed, threads=threads)
        snakemakeutils.dump_object(out.fq_out.to_fq_dict(), Path(output.FASTQ))
        with open(output.INFORMS, 'w') as handle:
            json.dump(out.informs, handle, indent=2)

rule preprocess_map_reads_post:
    """
    Maps the downsampled reads against the reference genome.
    """
    input:
        FASTA = iterativemapping.INPUT_FASTA_REF,
        IO_FASTQ = rules.preprocess_downsample_by_segment.output.FASTQ
    output:
        BAM = 'preprocess/mapping_post/bam.io',
        JSON = 'preprocess/mapping_post/stats.json'
    params:
        dir_ = lambda wildcards: 'preprocess/mapping_post',
        input_type = config['input_type']
    threads: 8
    run:
        from camel.app.components.workflows.utils.fastqinput import FastqInput
        from camel.app.io.tooliofile import ToolIOFile
        from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow

        # Run the workflow
        workflow = ReadMappingWorkflow(Path(str(params.dir_)).absolute())
        out = workflow.run(
            fastq_in=FastqInput.from_fq_dict(Path(input.IO_FASTQ), params.input_type),
            fasta_ref=snakemakeutils.load_object(Path(str(input.FASTA)))[0].path,
            threads=threads
        )

        # Save output
        snakemakeutils.dump_object([ToolIOFile(out.path_bam)], Path(output.BAM))
        with open(output.JSON, 'w') as handle:
            json.dump(out.stats, handle, indent=2)

rule preprocess_stats:
    """
    Combines the stats for the downsampling by segment.
    """
    input:
        JSON_pre = rules.preprocess_map_reads_pre.output.JSON,
        JSON_post = rules.preprocess_map_reads_post.output.JSON
    output:
        TSV = 'preprocess/stats.tsv'
    run:
        import pandas as pd
        with open(input.JSON_pre) as handle:
            data_pre = json.load(handle)
        with open(input.JSON_post) as handle:
            data_post = json.load(handle)

        records_out = []
        for seq_id in data_pre['by_chr']:
            records_out.append({
                'segment': seq_id.split('-')[-1],
                'depth_median_pre': data_pre['by_chr'][seq_id]['depth_median'],
                'depth_median_post': data_post['by_chr'][seq_id]['depth_median'],
                'covered_rate_pre': data_pre['by_chr'][seq_id]['covered_rate'],
                'covered_rate_post': data_post['by_chr'][seq_id]['covered_rate'],
            })
        data_out = pd.DataFrame(records_out)
        data_out.sort_values(by='segment', inplace=True)
        data_out.to_csv(output.TSV, sep='\t', index=False)

rule preprocess_report:
    """
    Creates a report for the pre-processing workflow.
    """
    input:
        TSV = rules.preprocess_stats.output.TSV,
        FASTA = iterativemapping.INPUT_FASTA_REF,
        BAM = rules.preprocess_map_reads_post.output.BAM
    output:
        VAL_HTML = 'preprocess/report/html.iob' # preprocess.OUTPUT_REPORT
    params:
        dir_ = 'preprocess/report',
        name = config['sample_name'],
        max_depth = config['downsampling'].get('coverage_max_by_segment', 250),
        gap_depth_cutoff = config['low_depth'].get('gap_depth_cutoff', 5)
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.viral_consensus.reportersegmentdownsampling import ReporterSegmentDownsampling
        reporter = ReporterSegmentDownsampling()
        reporter.add_input_files({'TSV': [ToolIOFile(Path(input.TSV))]})
        snakemakeutils.add_pickle_inputs(reporter, input, keys=['BAM', 'FASTA'])
        reporter.update_parameters(
            max_depth=params.max_depth, gap_depth_cutoff=params.gap_depth_cutoff, name=str(params.name))
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

#####################
# Amplicon clipping #
#####################
rule preprocess_report_ampliconclip:
    """
    Creates the report for the amplicon clipping.
    """
    input:
        INFORMS = rules.preprocess_downsample_by_segment.output.INFORMS
    output:
        HTML = 'preprocess/ampliconclip/html.iob' # preprocess.OUTPUT_CLIPPING_REPORT
    params:
        dir_ = 'preprocess/ampliconclip'
    run:
        from camel.app.tools.pipelines.viral_consensus.reporterampliconclip import ReporterAmpliconClip
        reporter = ReporterAmpliconClip()
        with open(input.INFORMS) as handle:
            informs_amplicon_clip = next(d for d in json.load(handle) if 'ampliconclip' in d['_name'])
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        reporter.add_input_informs({'ampliconclip': informs_amplicon_clip})
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule preprocess_report_ampliconclip_empty:
    """
    Creates and empty report for the amplicon clipping.
    """
    output:
        HTML = 'preprocess/ampliconclip/html-empty.iob' # preprocess.OUTPUT_CLIPPING_REPORT_EMPTY
    params:
        dir_ = 'preprocess/ampliconclip'
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Primer removal', Path(output.HTML))

#########
# Other #
#########
rule preprocess_combine_informs:
    """
    Combines the informs for the pre-processing workflow.
    """
    input:
        INFORMS_ampligone = rules.preprocess_ampligone_to_bed.output.INFORMS if 'ampligone' in config['analyses'] else [],
        INFORMS_mapping = rules.preprocess_map_reads_pre.output.INFORMS,
        INFORMS_downsample = rules.preprocess_downsample_by_segment.output.INFORMS
    output:
        INFORMS = 'preprocess/report/informs.io' # preprocess.OUTPUT_INFORMS
    run:
        data_out = []
        # Ampligone informs (if executed)
        if len(input.INFORMS_ampligone) > 0:
            informs_ampligone = snakemakeutils.load_object(Path(input.INFORMS_ampligone))
            if isinstance(informs_ampligone, dict):
                data_out.append(informs_ampligone)
            else:
                data_out.append(informs_ampligone[0])
        # noinspection PyUnresolvedReferences
        for path_json in [Path(x) for x in (input.INFORMS_mapping, input.INFORMS_downsample)]:
            with path_json.open() as handle:
                informs_all = json.load(handle)
                for inform in informs_all:
                    inform['_tag'] = 'Pre-processing'
                    data_out.append(inform)
        snakemakeutils.dump_object(data_out, Path(output.INFORMS))

rule preprocess_create_summary_out:
    """
    Creates the summary output for the pre-processing workflow.
    """
    input:
        TSV = rules.preprocess_stats.output.TSV
    output:
        FILE = 'preprocess/report/summary_preprocess.{ext}' # preprocess.OUTPUT_PRE_PROCESS_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import pandas as pd

        # Extract metrics
        data_in = pd.read_table(input.TSV, na_values=['-'], keep_default_na=False)
        rows_out = []
        for row in data_in.to_dict('records'):
            rows_out.extend([
                [f"preprocess_{row['segment']}_{k}", row[k]]
            for k in ('depth_median_pre', 'depth_median_post', 'covered_rate_pre', 'covered_rate_post')])

        # Save to summary
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'preprocess')
