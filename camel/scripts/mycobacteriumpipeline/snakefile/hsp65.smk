from pathlib import Path

from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake.step import Step
from camel.snakefiles import gene_detection

rule hps65_report:
    """
    Creates a HTML report for the HSP65 (sub)species identification.
    """
    input:
        INFORMS_hits = str(gene_detection.OUTPUT_ALL_HITS).format(db='hsp65'),
        INFORMS_columns = str(gene_detection.OUTPUT_COLUMNS).format(db='hsp65'),
        FASTA_DB = str(gene_detection.GENE_DETECTION_FASTA).format(db='hsp65')
    output:
        VAL_HTML = 'hsp65/report/html.iob', # hsp65.OUTPUT_REPORT
        INFORMS = 'hsp65/report/informs.io'
    params:
        dir_ = 'hsp65/report',
        hit_type = config['gene_detection']['options']['method']
    run:
        from camel.app.tools.pipelines.mycobacterium.hsp65reporter import Hsp65Reporter
        reporter = Hsp65Reporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        reporter.update_parameters(hit_type=params.hit_type)
        snakemakeutils.add_io_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule hps65_report_empty:
    """
    Creates an empty HTML report for the HSP65 (sub)species identification.
    """
    output:
        VAL_HTML = 'hsp65/report/html-empty.iob' # hsp65.OUTPUT_REPORT_EMPTY
    params:
        dir_ = 'hsp65/report'
    run:
        from camelcore.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.hsp65reporter import Hsp65Reporter
        section = Hsp65Reporter.generate_empty_section()
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule hsp65_dump_summary_info:
    """
    Dumps the summary information for the hsp65 workflow in tabular format.
    """
    input:
        INFORMS = rules.hps65_report.output.INFORMS
    output:
        FILE = 'hsp65/summary/summary.{ext}' # hsp65.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [('hsp65_species', informs['hits']),]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'hsp65')
