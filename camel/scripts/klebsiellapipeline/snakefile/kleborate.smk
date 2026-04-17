from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly

rule kleborate_run:
    """
    Runs the Kleborate tool.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TSV = 'kleborate/tool/tsv.io',
        INFORMS = 'kleborate/tool/informs.io'
    params:
        dir_ = 'kleborate/tool'
    run:
        from camel.app.tools.pipelines.klebsiella.kleborate import Kleborate
        kleborate = Kleborate()
        kleborate.update_parameters(all=True)
        snakemakeutils.add_io_inputs(kleborate, input)
        step = Step(rule_name=str(rule), tool=kleborate, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(kleborate, output)

rule kleborate_reporter:
    """
    Creates the output report for Kleborate.
    """
    input:
        TSV = rules.kleborate_run.output.TSV,
        INFORMS_kleborate = rules.kleborate_run.output.INFORMS
    output:
        HTML = 'kleborate/report/html.iob' # kleborate.OUTPUT_REPORT
    params:
        dir_ = 'kleborate/report'
    run:
        from camel.app.tools.pipelines.klebsiella.kleboratereporter import KleborateReporter
        reporter = KleborateReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule kleborate_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'kleborate/report/html-empty.iob' # kleborate.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Kleborate', Path(output.VAL_HTML))

rule kleborate_create_summary:
    """
    Creates a tabular summary output for Kleborate.
    """
    input:
        INFORMS = rules.kleborate_run.output.INFORMS
    output:
        FILE = 'kleborate/summary/summary_kleborate.{ext}'
    params:
        keys_kept = [
            'ST', 'Yersiniabactin', 'YbST', 'wzi', 'K_locus', 'K_locus_identity', 'virulence_score',
            'resistance_score', 'O_locus', 'O_locus_identity'],
        ext = lambda wildcards: wildcards.ext
    run:
        informs_in = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [(f'kleborate_{key}', informs_in[key]) for key in params.keys_kept]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'kleborate')
