from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection


rule sccmec_typing_run:
    """
    Determines the SCCmec information based on the detected BLAST hits.
    """
    input:
        YAML = config['sccmec_typing']['profiles'],
        VAL_HITS = str(gene_detection.OUTPUT_ALL_HITS).format(db='sccmec_genes')
    output:
        INFORMS = 'sccmec_typing/tool/informs.io'
    params:
        dir_ = 'sccmec_typing/tool'
    run:
        from camel.app.tools.sccmectyping.sccmectyping import SCCmecTyping
        from camel.app.core.io.tooliofile import ToolIOFile
        sccmec_typing = SCCmecTyping()
        sccmec_typing.add_input_files({'YML': [ToolIOFile(Path(input.YAML))]})
        snakemakeutils.add_pickle_input(sccmec_typing, 'VAL_HITS', Path(input.VAL_HITS))
        step = Step(rule_name=str(rule), tool=sccmec_typing, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(sccmec_typing, output)

rule sccmec_typing_report:
    """
    Creates a report section for the SCCmec typing.
    """
    input:
        INFORMS = rules.sccmec_typing_run.output.INFORMS
    output:
        HTML = 'sccmec_typing/report/html.iob' # sccmectyping.OUTPUT_REPORT
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.core.reports.htmlreportsection import HtmlReportSection
        section = HtmlReportSection('SCC<i>mec</i> type', 3)
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        section.add_table(
            [[f"{complex_['name']}:", complex_['value'] if complex_['value'] is not None else '-'] for
             complex_ in informs['complexes']], table_attributes=[('class', 'information')])
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule sccmec_typing_report_empty:
    """
    Creates an empty report when SCCmec typing is disabled.
    """
    output:
        VAL_HTML = 'sccmec_typing/report/html-empty.iob' # sccmectyping.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('SCC<i>mec</i> type', Path(output.VAL_HTML))

rule sccmec_typing_summary:
    """
    Creates the summary output for the SCCmec workflow.
    """
    input:
        INFORMS = rules.sccmec_typing_run.output.INFORMS
    output:
        FILE = 'sccmec_typing/summary_sccmec.{ext}' # sccmectyping.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [
            (complex_['key'], complex_['value'] if complex_['value'] is not None else '-')
            for complex_ in informs['complexes']
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'sccmec_typing')
