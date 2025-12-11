import shutil
from pathlib import Path

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import human_read_scrubbing
from camel.snakefiles.human_read_scrubbing import get_removed


rule scrubbing_fasta_fa2fq:
    """
    Convert the input FASTA file to FASTQ format to ensure compatibility with the scrubber tool, which accepts only 
    FASTQ input.
    """
    input:
        FASTA = human_read_scrubbing.INPUT_FASTA
    output:
        FASTQ = 'human_read_scrubbing/fasta/input/fastq.io'
    params:
        dir_ = 'human_read_scrubbing/fasta/input'
    run:
        from camel.app.core.utils import fastautils
        fasta_path_in = (snakemakeutils.load_object(Path(input.FASTA)))[0].path
        fastq_path_out = Path(str(params.dir_), f"{fasta_path_in.stem}.fastq")
        fastautils.convert_fasta_to_fastq(fasta_path_in, fastq_path_out)
        snakemakeutils.dump_object([ToolIOFile(fastq_path_out)], Path(output.FASTQ))

rule scrubbing_decompress_fastq_se:
    """
    Decompresses the input SE FASTQ file (if needed). 
    """
    input:
        FASTQ = str(human_read_scrubbing.INPUT_FASTQ).format(input_format='fastq_se')
    output:
        FASTQ = 'human_read_scrubbing/fastq_se/decompress/fastq.io'
    params:
        dir_ = 'human_read_scrubbing/fastq_se/decompress',
        name = config['input']['sample_name']
    threads: 4
    run:
        from camel.app.core.utils import fileutils
        path_fq = snakemakeutils.load_object(Path(input.FASTQ))[0].path
        if not fileutils.is_gzipped(path_fq):
            logger.info(f'Input FASTQ file is already decompressed')
            path_out = path_fq
        else:
            path_out = Path(params.dir_, f'{params.name}.fastq').absolute()
            fileutils.gzip_extract(path_fq, path_out, threads=threads)
        snakemakeutils.dump_object([ToolIOFile(path_out)], Path(output.FASTQ))

rule scrubbing_interleave_fastq_pe:
    """
    Interleaves the input PE FASTQ data.
    """
    input:
        FASTQ = lambda wildcards: human_read_scrubbing.INPUT_FASTQ.format(input_format='fastq_pe')
    output:
        FASTQ = 'human_read_scrubbing/fastq_pe/interleave_pe/fastq.io'
    params:
        dir_ = 'human_read_scrubbing/fastq_pe/interleave_pe'
    threads: 8
    run:
        from camel.app.tools.seqtk.seqtkmergepe import SeqtkMergePE
        merge_pe = SeqtkMergePE()
        merge_pe.add_input_files({
            'FASTQ_PE': snakemakeutils.load_object(Path(str(input.FASTQ)))
        })
        step = Step(rule_name=str(rule), tool=merge_pe, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(merge_pe, output)

rule scrubbing_select_input:
    """
    Selects the input for the read scrubbing tool.
    """
    input:
        FASTQ_PE_IL = lambda wildcards: rules.scrubbing_interleave_fastq_pe.output.FASTQ if wildcards.input_format == 'fastq_pe' else [],
        FASTQ_SE_GZ = lambda wildcards: rules.scrubbing_decompress_fastq_se.output.FASTQ if wildcards.input_format == 'fastq_se' else [],
        FASTQ_FORM_FASTA = lambda wildcards: rules.scrubbing_fasta_fa2fq.output.FASTQ if wildcards.input_format == 'fasta' else []
    output:
        FASTQ = 'human_read_scrubbing/{input_format}/scrubbing/in/fastq.io'
    params:
        input_format = lambda wildcards: wildcards.input_format
    run:
        if params.input_format == 'fastq_pe':
            shutil.copyfile(Path(str(input.FASTQ_PE_IL)), output.FASTQ)
        elif params.input_format == 'fastq_se':
            shutil.copyfile(Path(str(input.FASTQ_SE_GZ)), output.FASTQ)
        elif params.input_format == 'fasta':
            shutil.copyfile(Path(str(input.FASTQ_FORM_FASTA)), output.FASTQ)
        else:
            raise ValueError(f'Unsupported input format: {params.input_format}')

rule scrubbing_run_scrubber:
    """
    Runs the NCBI human read scrubber on the input FASTQ, only accepts single gunzipped FASTQ files.
    """
    input:
        FASTQ_SE = rules.scrubbing_select_input.output.FASTQ,
        DB = config['read_scrubbing']['db'] if config['read_scrubbing'].get('db') is not None else []
    output:
        FASTQ_SCRUBBED = 'human_read_scrubbing/{input_format}/scrubbing/fastq_scrubbed.io',
        FASTQ_REMOVED = 'human_read_scrubbing/{input_format}/scrubbing/fastq_removed.io',
        INFORMS = 'human_read_scrubbing/{input_format}/scrubbing/informs.io'
    params:
        running_dir = lambda wildcards: f'human_read_scrubbing/{wildcards.input_format}/scrubbing',
        input_format = lambda wildcards: wildcards.input_format,
        export_removed_reads = config.get('read_scrubbing', {}).get('export_removed_reads'),
        is_interleaved = lambda wildcards: True if wildcards.input_format == 'fastq_pe' else False,
        input_type = config['input']['type']
    threads: 8
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber

        scrubber = NcbiHumanReadScrubber()
        step = Step(rule_name=str(rule), tool=scrubber, dir_=Path(str(params.running_dir)))

        # Parameters
        scrubber.update_parameters(
            interleaved=bool(params.is_interleaved),
            outputfile='reads_kept.fastq',
            threads=threads)
        if params.export_removed_reads:
            outputfile_removed_reads = 'reads_removed.fastq'
            scrubber.update_parameters(
                export_human_reads=True,
                outputfile_removed=outputfile_removed_reads
            )

        # Add input files and run tool
        snakemakeutils.add_pickle_input(scrubber, 'FASTQ_SE', Path(input.FASTQ_SE))
        if len(input.DB) > 0:
            scrubber.add_input_files({'DB': [ToolIOFile(Path(input.DB))]})
        step.run()

        # Informs for hybrid input
        if params.input_type == 'hybrid':
            if params.input_format == 'fastq_pe':
                scrubber.informs['_tag'] = 'Illumina'
            else:
                scrubber.informs['_tag'] = 'ONT'

        # Store the output
        snakemakeutils.dump_tool_outputs(scrubber, output, keys=['FASTQ_SCRUBBED', 'INFORMS'])
        if 'FASTQ_REMOVED' not in scrubber.tool_outputs:
            snakemakeutils.dump_object([], Path(output.FASTQ_REMOVED))
        else:
            snakemakeutils.dump_object(scrubber.tool_outputs['FASTQ_REMOVED'], Path(output.FASTQ_REMOVED))

rule scrubbing_select_output_fastq:
    """
    Links the FASTQ SE output to the FASTQ output.
    """
    input:
        FASTQ = str(rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED).format(input_format='fastq_se')
    output:
        FASTQ = 'human_read_scrubbing/fastq_se/output/fastq.io'
    shell:
        """
        cp {input.FASTQ} {output.FASTQ}
        """

rule scrubbing_fasta_fq2fa:
    """
    Converts the FASTQ output from the scrubber to FASTA format. 
    """
    input:
        FASTQ_SCRUBBED = lambda wildcards: str(rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED).format(input_format='fasta')
    output:
        FASTA = 'human_read_scrubbing/fasta/output/fasta.io'
    params:
        dir_ = 'human_read_scrubbing/fasta/output'
    run:
        from Bio import SeqIO
        path_in = (snakemakeutils.load_object(Path(str(input.FASTQ_SCRUBBED))))[0].path
        path_out = Path(str(params.dir_), f"{path_in.stem}.fasta")

        with path_in.open('r') as fastq_file, path_out.open('w') as fasta_file:
            # Write the FASTQ SeqRecords in FASTA format, not using SeqIO.write directly because it writes multi FASTA
            fasta_out = SeqIO.FastaIO.FastaWriter(fasta_file, wrap=None)
            fasta_out.write_file(SeqIO.parse(fastq_file, 'fastq'))
        snakemakeutils.dump_object([ToolIOFile(path_out)], Path(output.FASTA))

rule scrubbing_fasta_fq2fa_removed:
    """
    Converts the FASTQ containing removed reads back to FASTA.
    """
    input:
        FASTQ_REMOVED = rules.scrubbing_run_scrubber.output.FASTQ_REMOVED,
        INFORMS_SCRUBBER = rules.scrubbing_run_scrubber.output.INFORMS
    output:
        FASTA_REMOVED = 'human_read_scrubbing/{input_format}/output/fasta_removed.io'
    params:
        running_dir = lambda wildcards: f'human_read_scrubbing/{wildcards.input_format}/output'
    run:
        from Bio import SeqIO

        scrubber_informs = snakemakeutils.load_object(Path(input.INFORMS_SCRUBBER))
        path_removed_in = (snakemakeutils.load_object(Path(input.FASTQ_REMOVED)))[0].path
        path_removed_out = Path(str(params.running_dir), f"{path_removed_in.stem}.fasta")

        if scrubber_informs['statistics']['count_removed'] != 0:
            with path_removed_in.open('r') as fastq_removed_file, path_removed_out.open('w') as fasta_removed_file:
                # Write the FASTQ SeqRecords in FASTA format, not using SeqIO.write directly because it writes multine fasta's
                fasta_removed_out = SeqIO.FastaIO.FastaWriter(fasta_removed_file, wrap=None)
                fasta_removed_out.write_file(SeqIO.parse(fastq_removed_file, 'fastq'))
            snakemakeutils.dump_object([ToolIOFile(path_removed_out)], Path(output.FASTA_REMOVED))
        else:
            snakemakeutils.dump_object([], Path(output.FASTA_REMOVED))

rule scrubbing_deinterleave_fastq_pe:
    """
    De-interleaves the PE FASTQ file after scrubbing.
    """
    input:
        FASTQ = lambda wildcards: {
            'scrubbed': str(rules.scrubbing_run_scrubber.output.FASTQ_SCRUBBED).format(input_format='fastq_pe'),
            'removed': str(rules.scrubbing_run_scrubber.output.FASTQ_REMOVED).format(input_format='fastq_pe'),
        }[wildcards.group]
    output:
        FASTQ = 'human_read_scrubbing/fastq_pe/output/{group}/fastq.io'
    params:
        dir_ = lambda wildcards: f'human_read_scrubbing/fastq_pe/output/{wildcards.group}',
        name = config['input']['sample_name']
    threads: 8
    run:
        from camel.app.tools.seqkit.seqkitsplit2 import SeqkitSplit2
        paths_fq = snakemakeutils.load_object(Path(str(input.FASTQ)))
        if len(paths_fq) == 0:
            snakemakeutils.dump_object([], Path(output.FASTQ))
        else:
            split2 = SeqkitSplit2()
            snakemakeutils.add_pickle_inputs(split2, input)
            split2.update_parameters(by_part=2)
            step = Step(rule_name=str(rule), tool=split2, dir_=Path(str(params.dir_)))
            step.run()
            snakemakeutils.dump_tool_outputs(split2, output)

rule scrubbing_fastq_gzip:
    """
    Compresses the FASTQ output.
    The 'group' wildcards corresponds to 'scrubbed' or 'removed' reads.
    Note: The FASTQ_PE reads need to be deinterleaved first.
    """
    input:
        FASTQ = lambda wildcards: \
            'human_read_scrubbing/{input_format}/scrubbing/fastq_{group}.io' if not wildcards.input_format == 'fastq_pe'
            else rules.scrubbing_deinterleave_fastq_pe.output.FASTQ
    output:
        FASTQ_GZ = 'human_read_scrubbing/{input_format}/compress/{group}/fastq_gz.io'
    threads: 4
    run:
        from camel.app.core.utils import fileutils
        output_io = []
        for io in snakemakeutils.load_object(Path(str(input.FASTQ))):
            path_out = io.path.parent / f'{io.path.name}.gz'
            fileutils.gzip_compress(io.path, path_out, threads=threads)
            output_io.append(ToolIOFile(path_out))
        snakemakeutils.dump_object(output_io, Path(output.FASTQ_GZ))

rule scrubbing_report:
    """
    Generates a tiny html report containing the tool name/version and a single phrase about how many reads were removed.
    """
    input:
        INFORMS_SCRUBBER = rules.scrubbing_run_scrubber.output.INFORMS,
        REMOVED = lambda wildcards: get_removed(input_format=wildcards.input_format) if config['read_scrubbing'].get('export_removed_reads') else []
    output:
        HTML = 'human_read_scrubbing/{input_format}/output/html.iob' # human_read_scrubbing.OUTPUT_REPORT
    params:
        running_dir = lambda wildcards: f'human_read_scrubbing/{wildcards.input_format}/output',
        input_format = lambda wildcards: wildcards.input_format,
        export_removed_reads = config.get('read_scrubbing', {}).get('export_removed_reads')
    run:
        from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubberreporter import NcbiHumanReadScrubberReporter

        reporter = NcbiHumanReadScrubberReporter()
        if params.export_removed_reads:
            snakemakeutils.add_pickle_inputs(reporter, input, keys=['REMOVED'])
        snakemakeutils.add_pickle_inputs(reporter, input, keys=['INFORMS_SCRUBBER'])
        reporter.update_parameters(input_format=str(params.input_format))
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.running_dir)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule scrubbing_create_summary:
    """
    Creates the tabular summary output for the HRRT assay.
    """
    input:
        INFORMS_SCRUBBER = rules.scrubbing_run_scrubber.output.INFORMS
    output:
        OUT = 'human_read_scrubbing/{input_format}/output/summary_out.{ext}' # human_read_scrubbing.OUTPUT_SUMMARY
    params:
        input_format = lambda wildcards: wildcards.input_format,
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS_SCRUBBER))

        subject_map = {
            'fastq_pe': 'read_pairs',
            'fastq_se': 'reads',
            'fasta': 'contigs'
        }
        subject = subject_map.get(str(params.input_format), 'contigs')
        count_total = informs['statistics']['count_total']
        count_removed = informs['statistics']['count_removed']

        data_summary = [
            ('scrubbing_tool_version', informs['_name_full']),
            (f'scrubbing_{subject}_in', count_total),
            (f'scrubbing_{subject}_removed', count_removed)
        ]
        snakemakeutils.export_summary(data_summary, Path(output.OUT), str(params.ext), 'human_read_scrubbing')

rule scrubbing_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'human_read_scrubbing/{input_format}/output/html-empty.iob' # human_read_scrubbing.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Human read removal', Path(output.VAL_HTML))
