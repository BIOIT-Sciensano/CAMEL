from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import lrefinder as lrefinder_workflow


rule run_lrefinder:
    """
    Runs the LRE-Finder tool.
    """
    input:
        IO = lrefinder_workflow.get_input(config)
    output:
        INFORMS = 'lrefinder/tool/informs.io' # lrefinder.OUTPUT_INFORMS
    params:
        dir_ = 'lrefinder/tool',
        input_type = config['input']['type']
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.lrefinder.lrefinder import LREFinder
        lrefinder = LREFinder()
        if params.input_type != 'fasta':
            key_reads = 'PE' if params.input_type == 'illumina' else 'SE'
            fq_dict = snakepipelineutils.extract_fq_input(Path(input.IO), key_se = 'FASTQ_SE', drop_empty = True, read_type = key_reads)
            lrefinder.add_input_files(fq_dict)
        else:
            # When the input type is FASTA the simulated reads are used as input
            snakemakeutils.add_pickle_input(lrefinder, 'FASTQ_PE', Path(input.IO))
        step = Step(rule_name=str(rule), tool=lrefinder, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(lrefinder, output)

rule lre_finder_report:
    """
    Creates the HTML output for the LRE-Finder tool.
    """
    input:
        INFORMS_lrefinder = rules.run_lrefinder.output.INFORMS
    output:
        HTML = 'lrefinder/report/html.iob' # lrefinder_workflow.OUTPUT_REPORT
    params:
        dir_ = 'lrefinder/report',
        input_type = config['input']['type']
    run:
        from camel.app.tools.lrefinder.lrefinderreporter import LREFinderReporter
        reporter = LREFinderReporter()
        snakemakeutils.add_pickle_inputs(reporter, input)
        if params.input_type == 'fasta':
            reporter.update_parameters(pseudo_reads=True)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule lre_finder_report_empty:
    """
    Creates an empty HTML output when the LRE-Finder tool is disabled.
    """
    output:
        HTML = 'lrefinder/report/html_empty.iob' # lrefinder_workflow.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('LRE-Finder', Path(output.HTML))

rule lre_finder_summary:
    """
    Creates the summary output for the LRE-Finder tool.
    """
    input:
        INFORMS = rules.run_lrefinder.output.INFORMS
    output:
        FILE = 'lrefinder/summary/summary.{ext}' # lrefinder_workflow.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [(f'lrefinder_{key}', informs[key]) for key in ('species', 'genes', 'mutations')]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'lrefinder')
