from pathlib import Path

from camel.app.error import PipelineExecutionError
from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import trimming_ont


rule trimming_ont_nanoplot_pre:
    """
    Creates NanoPlot reports for the raw reads. 
    """
    input:
        FASTQ = trimming_ont.INPUT_ONT_FASTQ
    output:
        HTML = 'trimming_ont/nanoplot-pre/html.io', # trimming_ont.OUTPUT_NANOPLOT_HTML_PRE,
        TSV = 'trimming_ont/nanoplot-pre/txt.io', # trimming_ont.OUTPUT_NANOPLOT_TXT_PRE,
        INFORMS = 'trimming_ont/nanoplot-pre/informs.io' # trimming_ont.OUTPUT_NANOPLOT_INFORMS_PRE
    params:
        dir_ = 'trimming_ont/nanoplot-pre'
    threads: 4
    run:
        from camel.app.tools.nanoplot.nanoplot import NanoPlot
        nanoplot = NanoPlot()
        snakemakeutils.add_pickle_inputs(nanoplot, input)
        step = Step(rule_name=str(rule), tool=nanoplot, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(nanoplot, output)

rule trimming_ont_seqkit:
    """
    Read trimming using seqkit.
    """
    input:
        FASTQ = trimming_ont.INPUT_ONT_FASTQ
    output:
        FASTQ = 'trimming_ont/seqkit/fastq.io', # trimming_ont.OUTPUT_READS,
        INFORMS = 'trimming_ont/seqkit/informs.io' # trimming_ont.OUTPUT_INFORMS
    threads: 4
    priority: 1
    params:
        dir_ = 'trimming_ont/seqkit',
        sample_name = config.get('sample_name', 'reads'),
        min_length = config.get('read_trimming', {}).get('ont', {}).get('min_length', 500),
        min_qual = config.get('read_trimming', {}).get('ont', {}).get('min_qual', 7)
    run:
        from camel.app.tools.seqkit.seqkitseq import SeqkitSeq
        seqkit = SeqkitSeq()
        seqkit.update_parameters(
            min_length=params.min_length,
            min_qual=params.min_qual,
            output_filename=f'{params.sample_name}-filtered.fastq',
            threads=threads
        )
        snakemakeutils.add_pickle_inputs(seqkit, input)
        step = Step(rule_name=str(rule), tool=seqkit, dir_=Path(params.dir_))
        step.run()
        if seqkit.informs['nb_seqs_out'] == 0:
            raise PipelineExecutionError('No reads left after filtering.')
        snakemakeutils.dump_tool_outputs(seqkit, output)

rule trimming_ont_nanoplot_post:
    """
    Creates NanoPlot reports of the trimmed reads.
    """
    input:
        FASTQ = rules.trimming_ont_seqkit.output.FASTQ
    output:
        HTML = 'trimming_ont/nanoplot-post/html.io', # trimming_ont.OUTPUT_NANOPLOT_HTML_POST,
        TSV = 'trimming_ont/nanoplot-post/txt.io', # trimming_ont.OUTPUT_NANOPLOT_TXT_POST,
        INFORMS = 'trimming_ont/nanoplot-post/informs.io' # trimming_ont.OUTPUT_NANOPLOT_INFORMS_POST
    params:
        dir_ = 'trimming_ont/nanoplot-post'
    threads: 4
    run:
        from camel.app.tools.nanoplot.nanoplot import NanoPlot
        nanoplot = NanoPlot()
        snakemakeutils.add_pickle_inputs(nanoplot, input)
        step = Step(rule_name=str(rule), tool=nanoplot, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(nanoplot, output)

rule trimming_ont_report:
    """
    Creates the HTML report with the trimming output files and statistics.
    """
    input:
        HTML_PRE = rules.trimming_ont_nanoplot_pre.output.HTML,
        HTML_POST = rules.trimming_ont_nanoplot_post.output.HTML,
        FASTQ_PE = rules.trimming_ont_seqkit.output.FASTQ,
        INFORMS_trimming = rules.trimming_ont_seqkit.output.INFORMS,
        INFORMS_nanoplot_pre = rules.trimming_ont_nanoplot_pre.output.INFORMS,
        INFORMS_nanoplot_post= rules.trimming_ont_nanoplot_post.output.INFORMS
    output:
        VAL_HTML = 'trimming_ont/report/html.iob' # trimming_ont.OUTPUT_REPORT
    params:
        dir_ = 'trimming_ont/report',
        export_fastq = config['read_trimming'].get('export_fastq', 'false')
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimmingont import ReporterTrimmingONT
        reporter = ReporterTrimmingONT()
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        reporter.update_parameters(export_fastq=str(params.export_fastq))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule trimming_ont_dump_summary_info:
    """
    Dumps the summary information from the read trimming pipeline.
    """
    input:
        INFORMS_trimming = rules.trimming_ont_seqkit.output.INFORMS,
        INFORMS_nanoplot = rules.trimming_ont_nanoplot_pre.output.INFORMS
    output:
        FILE = 'trimming_ont/summary/summary_out.{ext}' # trimming_ont.OUTPUT_SUMMARY
    params:
        dir_ = 'trimming_ont/summary',
        ext = lambda wildcards: wildcards.ext
    run:
        informs_filtering = snakemakeutils.load_object(Path(input.INFORMS_trimming))
        informs_nanoplot = snakemakeutils.load_object(Path(input.INFORMS_nanoplot))
        summary_data = [
            ('trim_ont_reads_in', informs_filtering['nb_seqs_in']),
            ('trim_ont_reads_out', informs_filtering['nb_seqs_out']),
            ('trim_ont_min_length', informs_filtering['min_length']),
            ('trim_ont_min_qual', informs_filtering['min_qual']),
            ('trim_ont_median_read_length', informs_nanoplot['median_read_length']),
            ('trim_ont_mean_qual', informs_nanoplot['mean_qual']),
            ('trim_ont_median_qual', informs_nanoplot['median_qual']),
            ('trim_ont_tool_version', informs_filtering['_name'])
        ]
        if params.ext == 'json':
            # informs_duplicates = snakemakeutils.load_object(Path(input.INFORMS_duplicates))
            summary_data.append(('duplication_rate', 'TODO!'))
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), 'trimming_ont')

rule trimming_ont_to_dict:
    """
    Combines the trimmed reads into a dictionary.
    """
    input:
        FASTQ = rules.trimming_ont_seqkit.output.FASTQ
    output:
        IO = 'trimming_ont/fastq_all.io' # trimming_ont.OUTPUT_DICT
    run:
        output_dict = {
            'SE': snakemakeutils.load_object(Path(input.FASTQ))
        }
        snakemakeutils.dump_object(output_dict,Path(output.IO))
