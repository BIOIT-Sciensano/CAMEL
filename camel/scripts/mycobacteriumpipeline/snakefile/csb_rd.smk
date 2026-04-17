from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection


rule rd_csb_report:
    """
    This rule generates the output HTML report for the RD - csb species identification.
    """
    input:
        FASTA = str(gene_detection.GENE_DETECTION_FASTA).format(db='csb_rd'),
        HITS = str(gene_detection.OUTPUT_ALL_HITS).format(db='csb_rd'),
        INFORMS_columns = str(gene_detection.OUTPUT_COLUMNS).format(db='csb_rd')
    output:
        VAL_HTML = 'csb_rd/report/html.iob', # csb_rd.OUTPUT_CSB_RD_REPORT
        INFORMS = 'csb_rd/report/informs.io'
    params:
        dir_ = 'csb_rd/report',
        detection_method = config['gene_detection']['options']['method'],
        input_type = config['input']['type']
    run:
        from camel.app.tools.pipelines.mycobacterium.rdcsbreporter import RdCsbReporter
        reporter = RdCsbReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        reporter.update_parameters(hit_type=params.detection_method)
        if params.input_type in ('fasta', 'fasta_with_vcf'):
            reporter.update_parameters(pseudo_reads=True)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule rd_csb_report_empty:
    """
    This generates an empty RD1 / csb report. 
    """
    output:
        VAL_HTML = 'csb_rd/report/html-empty.iob' # csb_rd.OUTPUT_CSB_RD_REPORT_EMPTY
    params:
        dir_ = 'csb_rd/report'
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.mycobacterium.rdcsbreporter import RdCsbReporter
        section = RdCsbReporter.generate_empty_section()
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule rd_csb_dump_summary_info:
    """
    Dumps the summary information from the csb / RD1 workflow in tabular format.
    """
    input:
        INFORMS_report = rules.rd_csb_report.output.INFORMS
    output:
        FILE = 'csb_rd/summary.{ext}' # csb_rd.OUTPUT_CSB_RD_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS_report))
        data_summary = [
            ('csb_detected', 'csb' in informs['loci_detected']),
            ('RD1_detected', 'RD1' in informs['loci_detected']),
            ('RD9_detected', 'RD9' in informs['loci_detected']),
            ('rd_csb_species', informs['species'])
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'csb_rd')
