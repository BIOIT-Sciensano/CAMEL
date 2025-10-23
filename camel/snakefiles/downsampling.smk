from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import downsampling


rule downsampling_read_stats:
    """
    Determines the number of reads and bases in the input FASTQ files.
    """
    input:
        FASTQ = downsampling.INPUT_FASTQ
    output:
        INFORMS = 'downsampling/{read_key}/read_stats/informs.io'
    params:
        dir_ = lambda wildcards: f'downsampling/{wildcards.read_key}/read_stats'
    run:
        from camel.app.tools.seqtk.seqtksize import SeqtkSize
        seqtk_size = SeqtkSize()
        snakemakeutils.add_pickle_inputs(seqtk_size, input)
        step = Step(rule_name=str(rule), tool=seqtk_size, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_object(seqtk_size.informs, Path(output.INFORMS))

rule downsampling_calculate:
    """
    Calculates if downsampling is required and associated statistics.
    """
    input:
        INFORMS_stats = rules.downsampling_read_stats.output.INFORMS
    output:
        JSON = 'downsampling/{read_key}/calculate_stats/json.io',
        INFORMS = 'downsampling/{read_key}/calculate_stats/informs.io'
    params:
        dir_ = lambda wildcards: f'downsampling/{wildcards.read_key}/calculate_stats',
        fasta_ref = config['reference'].get('fasta'),
        expected_size = config['reference'].get('expected_size'),
        cov_target = config['downsampling'].get('coverage_max'),
        is_disabled = config['downsampling'].get('disabled', False),
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        from camel.app.core.utils import fastautils
        from camel.app.tools.pipelines.downsampling.downsamplecalculation import DownsampleCalculation
        ds_calc = DownsampleCalculation()
        step = Step(rule_name=str(rule), tool=ds_calc, dir_=Path(str(params.dir_)))

        # Determine the expected size
        if params.is_disabled is True:
            size_ref = None
        elif params.expected_size is not None:
            size_ref = int(params.expected_size)
        elif params.fasta_ref is not None:
            size_ref = fastautils.count_bases(params.fasta_ref)
        else:
            raise ValueError(
                "Unable to determine the expected size. Please specify either 'expected_size' or 'fasta' in the "
                "'reference' section of the config file.")

        snakemakeutils.add_pickle_inputs(ds_calc, input)
        ds_calc.update_parameters(
            size_ref_genome=size_ref, cov_target=params.cov_target, is_paired=bool(params.is_paired))
        step.run()
        snakemakeutils.dump_tool_outputs(ds_calc, output)

rule downsampling_seqtk:
    """
    Performs downsampling with seqtk (if needed).
    """
    input:
        FASTQ = downsampling.INPUT_FASTQ,
        INFORMS = rules.downsampling_calculate.output.INFORMS
    output:
        FASTQ = 'downsampling/{read_key}/seqtk/fastq.io',
        INFORMS = 'downsampling/{read_key}/seqtk/informs.io'
    params:
        dir_ = lambda wildcards: f"downsampling/{wildcards.read_key}/seqtk",
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        import logging
        import shutil
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample

        # Create working directory
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)

        # Parse statistics
        data_ds = snakemakeutils.load_object(Path(input.INFORMS))['stats']
        if data_ds['downsample_factor'] is None:
            # No downsampling -> create empty output files
            logging.info("No downsampling required, skipping seqtk")
            shutil.copyfile(input.FASTQ, output.FASTQ)
            snakemakeutils.dump_object([], Path(output.INFORMS))
        else:
            # Downsampling -> run seqtk
            seqtk = SeqtkSubsample()
            key_fastq = 'FASTQ_PE' if params.is_paired else 'FASTQ'
            seqtk.add_input_files({key_fastq: snakemakeutils.load_object(Path(input.FASTQ))})
            step = Step(rule_name=str(rule), tool=seqtk, dir_=dir_)
            seqtk.update_parameters(fraction=float(f"{data_ds['downsample_factor']:.6f}"))
            step.run()
            snakemakeutils.dump_object(seqtk.tool_outputs[key_fastq], Path(output.FASTQ))
            snakemakeutils.dump_object(seqtk.informs, Path(output.INFORMS))

rule downsampling_report:
    """
    Creates the downsampling report.
    """
    input:
        INFORMS_stats = rules.downsampling_calculate.output.INFORMS,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        HTML = 'downsampling/{read_key}/report/html.iob' # downsampling.OUTPUT_REPORT
    params:
        dir_ = lambda wildcards: f"downsampling/{wildcards.read_key}/report",
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe'
    run:
        from camel.app.tools.pipelines.downsampling.reporterdownsampling import ReporterDownsampling
        reporter = ReporterDownsampling()
        reporter.add_input_informs({'stats': snakemakeutils.load_object(Path(input.INFORMS_stats))})
        reporter.update_parameters(is_paired=bool(params.is_paired))
        informs_seqtk = snakemakeutils.load_object(Path(input.INFORMS_seqtk))
        if len(informs_seqtk) > 0:
            reporter.add_input_informs({'seqtk': informs_seqtk})
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_object(reporter.tool_outputs['HTML'], Path(output.HTML))

rule downsampling_summary_out:
    """
    Exports the summary information for the downsampling workflow.
    """
    input:
        INFORMS_stats = rules.downsampling_calculate.output.INFORMS,
        INFORMS_seqtk = rules.downsampling_seqtk.output.INFORMS
    output:
        FILE = 'downsampling/{read_key}/summary/summary_downsampling.{ext}' # downsampling.OUTPUT_SUMMARY
    params:
        is_paired = lambda wildcards: wildcards.read_key == 'fastq_pe',
        read_key = lambda wildcards: wildcards.read_key,
        add_read_key_prefix = config['input_type'] == 'hybrid',
        ext = lambda wildcards: wildcards.ext
    run:
        # Parse calculation output
        data_ds = snakemakeutils.load_object(Path(input.INFORMS_stats))['stats']

        # Add the number of output reads / read pairs to the summary output
        key_nb_reads = 'nb_read_pairs' if params.is_paired else 'nb_reads'
        key_out = f'{key_nb_reads}_out'
        if data_ds['downsample_factor'] is None:
            data_ds[key_out] = data_ds[f'{key_nb_reads}_in']
            data_ds['downsample_factor'] = '-'
        else:
            informs_in = snakemakeutils.load_object(Path(input.INFORMS_seqtk))
            if params.is_paired:
                data_ds[key_out] = informs_in['reads_count'] // 2
            else:
                data_ds[key_out] = informs_in['reads_count']

        # Construct the summary data
        summary_data = []
        for k, v in sorted(data_ds.items()):
            # Add read key information when multiple FASTQ datasets were processed
            summary_data.append((
                f'downsampling_{k}' if not params.add_read_key_prefix else f'downsampling_{params.read_key}_{k}',
                f'{v:.2f}' if isinstance(v, float) else str(v)))

        # Create output
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), f'downsampling_{params.read_key}')
