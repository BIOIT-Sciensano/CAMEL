from pathlib import Path

from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import trimming_illumina
from camel.app.core.snakemake.step import Step


rule trimming_illumina_fastqc_pre:
    """
    Creates FastQC reports for the raw reads. 
    """
    input:
        FASTQ = trimming_illumina.INPUT_FASTQ
    output:
        HTML = 'trimming_illumina/fastqc-pre/html.io', # trimming_illumina.OUTPUT_FASTQC_HTML_PRE
        TXT = 'trimming_illumina/fastqc-pre/txt.io' # trimming_illumina.OUTPUT_FASTQC_TXT_PRE
    params:
        dir_ = 'trimming_illumina/fastqc-pre'
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC()
        snakemakeutils.add_io_inputs(fastqc, input)
        step = Step(rule_name=str(rule), tool=fastqc, dir_=Path(params.dir_))
        fastqc.update_parameters(threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(fastqc, output)

rule trimming_illumina_trimmomatic:
    """
    Read trimming using trimmomatic.
    """
    input:
        FASTQ_PE = trimming_illumina.INPUT_FASTQ
    output:
        FASTQ_PE = 'trimming_illumina/trimmomatic/fastq-pe.io',
        FASTQ_SE_FORWARD = 'trimming_illumina/trimmomatic/fastq-se-fwd.io',
        FASTQ_SE_REVERSE = 'trimming_illumina/trimmomatic/fastq-se-rev.io',
        INFORMS = 'trimming_illumina/trimmomatic/informs.io'
    threads: 4
    priority: 1
    params:
        dir_ = 'trimming_illumina/trimmomatic',
        adapter = config.get('read_trimming', {}).get('adapter'),
        sample_name = config.get('sample_name', 'reads')
    run:
        from camel.app.tools.trimmomatic.trimmomatic import Trimmomatic
        trimmomatic = Trimmomatic()
        snakemakeutils.add_io_inputs(trimmomatic, input)
        trimmomatic.update_parameters(baseout=f'{params.sample_name}-trimmed.fastq.gz')
        if params.adapter is not None:
            trimmomatic.update_parameters(illuminaclip_PE=f'$TRIMMOMATIC_ADAPTER_DIR/{params.adapter}-PE.fa:2:30:10')
        step = Step(rule_name=str(rule), tool=trimmomatic, dir_=Path(params.dir_))
        trimmomatic.update_parameters(threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(trimmomatic, output)

rule trimming_illumina_fastp:
    """
    Read trimming using fastp.
    """
    input:
        FASTQ = trimming_illumina.INPUT_FASTQ
    output:
        FASTQ_PE = 'trimming_illumina/fastp/fastq-pe.io',
        FASTQ_SE_FWD = 'trimming_illumina/fastp/fastq-se-fwd.io',
        FASTQ_SE_REV = 'trimming_illumina/fastp/fastq-se-rev.io',
        HTML = 'trimming_illumina/fastp/html.iob',
        INFORMS = 'trimming_illumina/fastp/informs.io'
    threads: 4
    priority: 1
    params:
        dir_ = 'trimming_illumina/fastp',
        sample_name = config.get('sample_name', 'reads')
    run:
        from camel.app.tools.fastp.fastp import Fastp
        fastp = Fastp()
        snakemakeutils.add_io_inputs(fastp, input)
        fastp.update_parameters(
            output_name=f'{params.sample_name}-trimmed',
            # Adapter trimming
            detect_adapter_for_pe=True,
            # Leading
            cut_front=True,
            cut_front_window_size=1,
            cut_front_mean_quality=10,
            # Trailing
            cut_tail=True,
            cut_tail_window_size=1,
            cut_tail_mean_quality=10,
            # Sliding window
            cut_right=True,
            cut_right_window_size=4,
            cut_right_mean_quality=20,
            # Minimum length
            length_required=40,
            # Threads
            threads=threads
        )
        step = Step(rule_name=str(rule), tool=fastp, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(fastp, output)

rule trimming_illumina_fastqc_post:
    """
    Creates FastQC reports of the trimmed reads.
    """
    input:
        FASTQ = trimming_illumina.select_fastq_output(config)
    output:
        HTML = 'trimming_illumina/fastqc-post/html.io', # trimming_illumina.OUTPUT_FASTQC_HTML_POST
        TXT = 'trimming_illumina/fastqc-post/txt.io' # trimming_illumina.OUTPUT_FASTQC_TXT_POST
    params:
        dir_ = 'trimming_illumina/fastqc-post'
    threads: 4
    run:
        from camel.app.tools.fastqc.fastqc import FastQC
        fastqc = FastQC()
        snakemakeutils.add_io_inputs(fastqc, input)
        step = Step(rule_name=str(rule), tool=fastqc, dir_=Path(params.dir_))
        fastqc.update_parameters(threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(fastqc, output)

rule trimming_illumina_report_trimmomatic:
    """
    Creates the HTML report with the trimming output files and statistics for Trimmomatic.
    """
    input:
        HTML_PRE = rules.trimming_illumina_fastqc_pre.output.HTML,
        HTML_POST = rules.trimming_illumina_fastqc_post.output.HTML,
        FASTQ_PE = rules.trimming_illumina_trimmomatic.output.FASTQ_PE,
        FASTQ_SE_FORWARD = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_FORWARD,
        FASTQ_SE_REVERSE = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_REVERSE,
        INFORMS_trimming = rules.trimming_illumina_trimmomatic.output.INFORMS
    output:
        VAL_HTML = 'trimming_illumina/report/trimmomatic/html.iob'
    params:
        dir_ = 'trimming_illumina/report/trimmomatic',
        export_fastq = config['read_trimming'].get('export_fastq', False)
    run:
        from camel.app.tools.pipelines.read_trimming.reportertrimming import ReporterTrimming
        reporter = ReporterTrimming()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        reporter.update_parameters(export_fastq=params.export_fastq)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule trimming_illumina_report_fastp:
    """
    Creates the HTML report with the trimming output files and statistics for fastp.
    """
    input:
        HTML = rules.trimming_illumina_fastp.output.HTML,
        HTML_pre = rules.trimming_illumina_fastqc_pre.output.HTML,
        HTML_post = rules.trimming_illumina_fastqc_post.output.HTML,
        FASTQ_PE = rules.trimming_illumina_fastp.output.FASTQ_PE,
        FASTQ_SE_FWD = rules.trimming_illumina_fastp.output.FASTQ_SE_FWD,
        FASTQ_SE_REV = rules.trimming_illumina_fastp.output.FASTQ_SE_REV,
        INFORMS_fastp = rules.trimming_illumina_fastp.output.INFORMS
    output:
        VAL_HTML = 'trimming_illumina/report/fastp/html.iob'
    params:
        dir_ = 'trimming_illumina/report/fastp',
        export_fastq = config['read_trimming'].get('export_fastq', 'false')
    run:
        from camel.app.tools.fastp.fastpreporter import FastpReporter
        reporter = FastpReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        reporter.update_parameters(export_fastq=str(params.export_fastq))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule trimming_illumina_report_select:
    """
    Selects the report output based on the selected trimming method.
    """
    input:
        HTML = f"trimming_illumina/report/{config['read_trimming'].get('method', 'trimmomatic')}/html.iob"
    output:
        HTML = 'trimming_illumina/report/html.iob'
    shell:
        """
        cp {input.HTML} {output.HTML};
        """

rule trimming_illumina_summary_trimmomatic:
    """
    Dumps the summary information for trimming with trimmomatic. 
    """
    input:
        INFORMS = rules.trimming_illumina_trimmomatic.output.INFORMS
    output:
        FILE = 'trimming_illumina/summary/trimmomatic/summary.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs_trimmomatic = snakemakeutils.load_object(Path(input.INFORMS))
        summary_data = [
            ('trim_ilmn_pairs_in', informs_trimmomatic['paired_reads_in']),
            ('trim_ilmn_pairs_out', informs_trimmomatic['paired_reads_out'].split(' ')[0]),
            ('trim_ilmn_fwd_only_surviving', informs_trimmomatic['forward_only_reads'].split(' ')[0]),
            ('trim_ilmn_rev_only_surviving', informs_trimmomatic['reverse_only_reads'].split(' ')[0]),
            ('trim_ilmn_pairs_both_dropped', informs_trimmomatic['reads_drop'].split(' ')[0]),
            ('trim_ilmn_tool_version', informs_trimmomatic['_name_full'])
        ]
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), 'trimming_illumina')

rule trimming_illumina_summary_fastp:
    """
    Dumps the summary information for trimming with fastp.
    """
    input:
        INFORMS = rules.trimming_illumina_fastp.output.INFORMS
    output:
        FILE = 'trimming_illumina/summary/fastp/summary.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs_fastp = snakemakeutils.load_object(Path(input.INFORMS))
        summary_data = [
            ('trim_ilmn_pairs_in', informs_fastp['summary']['before_filtering']['total_read_pairs']),
            ('trim_ilmn_pairs_out', informs_fastp['summary']['after_filtering']['total_read_pairs']),
            ('trim_ilmn_bases_in', informs_fastp['summary']['before_filtering']['total_bases']),
            ('trim_ilmn_bases_out', informs_fastp['summary']['after_filtering']['total_bases']),
            ('trim_ilmn_q30_rate_in', informs_fastp['summary']['before_filtering']['q30_rate']),
            ('trim_ilmn_q30_rate_out', informs_fastp['summary']['after_filtering']['q30_rate']),
            ('trim_ilmn_duplication_rate', informs_fastp['duplication']['rate']),
            ('trim_ilmn_tool_version', informs_fastp['_name_full']),
        ]
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), 'trimming_illumina')

rule trimming_illumina_summary_select:
    """
    Select the summary output depending on the selected trimming method.
    """
    input:
        FILE = lambda wildcards: f"trimming_illumina/summary/{config['read_trimming'].get('method', 'trimmomatic')}/summary.{wildcards.ext}"
    output:
        FILE = 'trimming_illumina/summary/summary_trim.{ext}'
    params:
        dir_ = 'trimming_illumina/summary'
    shell:
        """
        cp {input.FILE} {output.FILE};
        """

rule trimming_illumina_to_dict_trimmomatic:
    """
    Combines the reads trimmed by Trimmomatic into a dictionary.
    """
    input:
        FASTQ_PE = rules.trimming_illumina_trimmomatic.output.FASTQ_PE,
        FASTQ_SE_FWD = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_FORWARD,
        FASTQ_SE_REV = rules.trimming_illumina_trimmomatic.output.FASTQ_SE_REVERSE
    output:
        IO = 'trimming_illumina/trimmomatic/fq_dict.io'
    run:
        output_dict = {'PE': snakemakeutils.load_object(Path(input.FASTQ_PE))}
        se_fwd = snakemakeutils.load_object(Path(input.FASTQ_SE_FWD))
        if len(se_fwd) > 0:
            output_dict['SE_FWD'] = se_fwd
        se_rev = snakemakeutils.load_object(Path(input.FASTQ_SE_REV))
        if len(se_rev) > 0:
            output_dict['SE_REV'] = se_rev
        snakemakeutils.dump_object(output_dict, Path(output.IO))

rule trimming_illumina_to_dict_fastp:
    """
    Combines the reads trimmed by fastp into a dictionary.
    """
    input:
        FASTQ_PE = rules.trimming_illumina_fastp.output.FASTQ_PE,
        FASTQ_SE_FWD = rules.trimming_illumina_fastp.output.FASTQ_SE_FWD,
        FASTQ_SE_REV = rules.trimming_illumina_fastp.output.FASTQ_SE_REV
    output:
        IO = 'trimming_illumina/fastp/fq_dict.io'
    run:
        output_dict = {'PE': snakemakeutils.load_object(Path(input.FASTQ_PE))}
        se_fwd = snakemakeutils.load_object(Path(input.FASTQ_SE_FWD))
        if len(se_fwd) > 0:
            output_dict['SE_FWD'] = se_fwd
        se_rev = snakemakeutils.load_object(Path(input.FASTQ_SE_REV))
        if len(se_rev) > 0:
            output_dict['SE_REV'] = se_rev
        snakemakeutils.dump_object(output_dict, Path(output.IO))

rule trimming_illumina_to_dict_select:
    """
    Selects the FASTQ dictionary based on the trimming method.
    """
    input:
        IO = f"trimming_illumina/{config['read_trimming'].get('method', 'trimmomatic')}/fq_dict.io"
    output:
        IO = trimming_illumina.OUTPUT_DICT
    shell:
        """
        cp {input.IO} {output.IO};
        """

rule trimming_illumina_informs_select:
    """
    Selects the informs based on the selected trimming method.
    """
    input:
        INFORMS = f"trimming_illumina/{config['read_trimming'].get('method', 'trimmomatic')}/informs.io"
    output:
        INFORMS = trimming_illumina.OUTPUT_INFORMS
    shell:
        """
        cp {input.INFORMS} {output.INFORMS};
        """
