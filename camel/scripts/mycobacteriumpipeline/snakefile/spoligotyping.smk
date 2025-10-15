from pathlib import Path

from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly, variant_calling


rule spoligotyping_downsample:
    """
    Optional downsampling if the coverage is too high.
    """
    input:
        IO_FASTQ = 'fq_dict.io',
        INFORMS_coverage = variant_calling.OUTPUT_DEPTH_INFORMS
    output:
        FASTQ_PE = 'spoligotyping/downsampling/fastq-ds.io',
        INFORMS_spoligo_param = 'spoligotyping/downsampling/informs-param.io'
    params:
        dir_ = 'spoligotyping/downsampling',
        read_key = 'SE' if config['input_type'] == 'ont' else 'PE'
    run:
        from camel.app.loggers import logger
        from camel.app.tools.seqtk.seqtksubsample import SeqtkSubsample
        estimated_coverage = snakemakeutils.load_object(Path(input.INFORMS_coverage))['median_depth']
        downsample_factor = 50 / max(estimated_coverage, 1)
        logger.info(f"Downsampling factor for SpoTyping: {downsample_factor:.2f}")
        if downsample_factor < 1:
            seqtk = SeqtkSubsample()
            seqtk.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO_FASTQ)))
            step = Step(rule_name=str(rule), tool=seqtk, dir_=Path(str(params.dir_)))
            seqtk.update_parameters(fraction=downsample_factor)
            step.run()
            snakemakeutils.dump_object({
                'min_strict': 5,
                'min_relaxed': 5,
                'downsample_factor': '{:.2f}'.format(downsample_factor)
            }, Path(output.INFORMS_spoligo_param))
            snakemakeutils.dump_object(seqtk.tool_outputs['FASTQ_PE'], Path(output.FASTQ_PE))
        else:
            fq_dict = snakemakeutils.load_object(Path(input.IO_FASTQ))
            snakemakeutils.dump_object(fq_dict[params.read_key], Path(output.FASTQ_PE))
            snakemakeutils.dump_object({
                'min_strict': max(round(estimated_coverage / 10), 3),
                'min_relaxed': max(round(estimated_coverage / 10), 3),
                'downsample_factor': 'NA'
            }, Path(output.INFORMS_spoligo_param))

rule spoligotyping_spotyping:
    """
    Runs the SpoTyping tool.
    """
    input:
        FASTQ = rules.spoligotyping_downsample.output.FASTQ_PE if config['input_type'] == 'illumina' else [],
        FASTA = assembly.OUTPUT_FASTA if config['input_type'] != 'illumina' else [],
        INFORMS_spoligo_param = rules.spoligotyping_downsample.output.INFORMS_spoligo_param if config['input_type'] == 'illumina' else []
    output:
        VAL_type_binary = 'spoligotyping/spotyping/VAL_binary.io',
        VAL_type_octal = 'spoligotyping/spotyping/VAL_octal.io',
        LOG = 'spoligotyping/spotyping/log.io',
        INFORMS = 'spoligotyping/spotyping/informs.io' # spoligotyping.OUTPUT_INFORMS
    params:
        dir_ = 'spoligotyping/spotyping',
        key = 'FASTQ' if config['input_type'] == 'illumina' else 'FASTA'
    run:
        from camel.app.tools.spotyping.spotyping import SpoTyping
        spotyping = SpoTyping()
        if params.key == 'FASTQ':
            snakemakeutils.add_pickle_input(spotyping, params.key, Path(input.FASTQ))
            spotyping_params = snakemakeutils.load_object(Path(input.INFORMS_spoligo_param))
            spotyping.update_parameters(
                swift='off', min_strict=spotyping_params['min_strict'], min_relaxed=spotyping_params['min_relaxed'])
        else:
            snakemakeutils.add_pickle_input(spotyping, params.key, Path(input.FASTA))
            spotyping.update_parameters(swift='off', fasta=True)
        step = Step(rule_name=str(rule), tool=spotyping, dir_=Path(str(params.dir_)))
        step.run()
        spotyping.informs['_tag'] = 'Spoligotyping'
        snakemakeutils.dump_tool_outputs(spotyping, output)

rule spoligotyping_report:
    """
    Creates a report for the spoligotyping.
    """
    input:
        VAL_type_binary = rules.spoligotyping_spotyping.output.VAL_type_binary,
        VAL_type_octal = rules.spoligotyping_spotyping.output.VAL_type_octal,
        LOG = rules.spoligotyping_spotyping.output.LOG,
        INFORMS_spotyping = rules.spoligotyping_spotyping.output.INFORMS,
        INFORMS_spoligo_param = rules.spoligotyping_downsample.output.INFORMS_spoligo_param if config['input_type'] == 'illumina' else []
    output:
        VAL_HTML = 'spoligotyping/report/html.iob', # spoligotyping.OUTPUT_REPORT
        INFORMS = 'spoligotyping/report/informs-report.io'
    params:
        dir_ = 'spoligotyping/report',
        input_type = config['input_type']
    run:
        from camel.app.tools.spotyping.spotypingreporter import SpoTypingReporter
        reporter = SpoTypingReporter()
        # noinspection PyUnresolvedReferences
        keys = [k for k, path in input.items() if len(path) > 0]
        snakemakeutils.add_pickle_inputs(reporter, input, keys=keys)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run_step()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule spoligotyping_report_empty:
    """
    Creates an empty report for the spoligotyping.
    """
    output:
        VAL_HTML = 'spoligotyping/report/html-empty.iob' # spoligotyping.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.spotyping.spotypingreporter import SpoTypingReporter
        section = SpoTypingReporter.generate_empty_section()
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule spoligotyping_dump_summary_info:
    """
    Dumps the summary information from the spoligotyping assay.
    """
    input:
        VAL_type_binary = rules.spoligotyping_spotyping.output.VAL_type_binary,
        VAL_type_octal = rules.spoligotyping_spotyping.output.VAL_type_octal,
        INFORMS = rules.spoligotyping_spotyping.output.INFORMS
    output:
        FILE = 'spoligotyping/summary/summary_out.{ext}' # spoligotyping.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [
            ('spoligotype_binary', snakemakeutils.load_object(Path(input.VAL_type_binary))[0].value),
            ('spoligotype_octal', snakemakeutils.load_object(Path(input.VAL_type_octal))[0].value),
            ('spoligotyping_tool_version', informs['_name_full']),
            ('sit_number', informs.get('metadata', {}).get('SIT', 'NA'))
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'spoligotyping')
