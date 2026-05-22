from pathlib import Path

from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection


rule serotype_detection_run:
    """
    Retrieves the serotype based on the H and O typing. 
    """
    input:
        HITS_O = str(gene_detection.OUTPUT_ALL_HITS).format(db='serotype_o'),
        HITS_H = str(gene_detection.OUTPUT_ALL_HITS).format(db='serotype_h')
    output:
        VAL_serotype = 'serotype_detection/tool/val-sero.io' # serotype_detection.OUTPUT_VAL
    params:
        dir_ =  'serotype_detection/tool'
    run:
        from camel.app.tools.pipelines.stec.serotypedetector import SerotypeDetectorEcoli
        from camel.app.core.snakemake.step import Step
        detector = SerotypeDetectorEcoli()
        snakemakeutils.add_io_inputs(detector, input)
        step = Step(rule_name=str(rule), tool=detector, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(detector, output)

rule serotype_detection_report:
    """
    Creates a simple report section for the detected serotype.
    """
    input:
        VAL_serotype = rules.serotype_detection_run.output.VAL_serotype if 'serotype' in config['analyses_selected'] else []
    output:
        HTML = 'serotype_detection/report/html.iob' # serotype_detection.OUTPUT_REPORT
    run:
        from camelcore.app.reports.htmlreportsection import HtmlReportSection
        from camelcore.app.io.tooliovalue import ToolIOValue
        if len(input['VAL_serotype']) > 0:
            serotype = snakemakeutils.load_object(Path(input.VAL_serotype))[0].value
        else:
            serotype = 'NA'
        section = HtmlReportSection(None)
        section.add_paragraph("Detected serotype: <b>{}</b>".format(serotype))
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule serotype_detection_dump_summary_info:
    """
    Dumps the summary information from the serotype detection.
    """
    input:
        VAL_serotype = rules.serotype_detection_run.output.VAL_serotype if 'serotype' in config['analyses_selected'] else []
    output:
         FILE = 'serotype_detection/summary/summary_out.{ext}' # serotype_detection.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        serotype = snakemakeutils.load_object(Path(input.VAL_serotype))[0].value
        data_summary = [('serotype', str(serotype))]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'serotype_determination')
