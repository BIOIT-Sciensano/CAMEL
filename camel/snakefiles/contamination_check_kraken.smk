from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import contamination_check_kraken


rule contamination_check_kraken2_run:
    """
    Assigns taxonomic labels to reads using KRAKEN2.
    """
    input:
        IO = str(contamination_check_kraken.get_input(config)),
        DB = config['contamination_check']['db']
    output:
        TSV = 'contamination_check/{input_format}/kraken2/tsv.io',
        TSV_report = 'contamination_check/{input_format}/kraken2/tsv-report.io',
        INFORMS = 'contamination_check/{input_format}/kraken2/informs.io'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/kraken2',
        input_format = lambda wildcards: wildcards.input_format
    threads: 8
    priority: 1
    run:
        from camel.app.core.io.tooliodirectory import ToolIODirectory
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.kraken.kraken2 import Kraken2

        kraken2 = Kraken2()
        if params.input_format == 'fastq_pe':
            kraken2.add_input_files(snakepipelineutils.extract_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        elif params.input_format == 'fastq_se':
            kraken2.add_input_files(snakepipelineutils.extract_fq_input(Path(input.IO), key_se='FASTQ', read_type='SE'))
        else:
            snakemakeutils.add_io_input(kraken2,'FASTA', Path(input.IO))
        kraken2.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        step = Step(rule_name=str(rule), tool=kraken2, dir_=Path(str(params.dir_)))
        kraken2.update_parameters(threads=threads)
        step.run()
        snakemakeutils.dump_io_outputs(kraken2, output)

rule contamination_check_kraken_report_parser:
    """
    Parses the Kraken report and looks for contamination at the species level. 
    """
    input:
        TSV = rules.contamination_check_kraken2_run.output.TSV_report,
        TSV_full = rules.contamination_check_kraken2_run.output.TSV
    output:
        TSV = 'contamination_check/{input_format}/kraken2/tsv-normalized.io',
        INFORMS = 'contamination_check/{input_format}/kraken2/informs-contamination.io'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/kraken2',
        expected_species = config['contamination_check']['expected_species'],
        allowed_species = config['contamination_check'].get('allowed_species', None),
        level_of_depth = config['contamination_check'].get('level_of_depth', 'S'),
        input_format = lambda wildcards: wildcards.input_format
    run:
        from camel.app.tools.kraken.krakenreportparser import KrakenReportParser
        report_parser = KrakenReportParser()
        snakemakeutils.add_io_inputs(report_parser, input)
        step = Step(rule_name=str(rule), tool=report_parser, dir_=Path(str(params.dir_)))
        report_parser.update_parameters(
            expected_species=params.expected_species if params.expected_species is not None else 'n/a')
        snakemakeutils.update_param_if_not_none(report_parser, 'level_of_depth', params)
        snakemakeutils.update_param_if_not_none(
            report_parser, 'allowed_species', params, transform=lambda x: ','.join(x))
        if params.input_format == 'fasta':
            report_parser.update_parameters(normalize_for_len=True)
        step.run()
        snakemakeutils.dump_io_outputs(report_parser, output)

rule contamination_check_krona:
    """
    Creates an interactive pie chart displaying the Kraken output.
    """
    input:
        TSV = rules.contamination_check_kraken2_run.output.TSV,
        DB = config['contamination_check'].get('db', [])
    output:
        HTML = 'contamination_check/{input_format}/krona/html.iob'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/krona'
    run:
        from camel.app.core.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.krona.krona import Krona
        krona = Krona()
        snakemakeutils.add_io_input(krona,'TSV', Path(input.TSV))
        if len(input.DB) > 0:
            krona.add_input_files({'DB': [ToolIODirectory(Path(input.DB, 'krona'))]})
        step = Step(rule_name=str(rule), tool=krona, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(krona, output)

rule contamination_check_report:
    """
    Creates a report containing the results of the contamination check.
    """
    input:
        HTML_Krona = rules.contamination_check_krona.output.HTML,
        INFORMS_species = rules.contamination_check_kraken_report_parser.output.INFORMS,
        INFORMS_kraken2 = rules.contamination_check_kraken2_run.output.INFORMS,
        TSV = rules.contamination_check_kraken_report_parser.output.TSV
    output:
        VAL_HTML = 'contamination_check/{input_format}/report/html.iob'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/report',
        input_format = lambda wildcards: wildcards.input_format,
        input_type = config['input']['type']
    run:
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination
        reporter = HtmlReporterContamination()
        snakemakeutils.add_io_inputs(reporter, input)

        # Update report title and output files for hybrid data
        if params.input_type == 'hybrid':
            reporter.update_parameters(
                suffix='illumina' if params.input_format == 'fastq_pe' else 'ont',
                suffix_title='Illumina' if params.input_format == 'fastq_pe' else 'ONT',
            )

        # Run the tool
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule contamination_check_report_empty:
    """
    Generates an empty contamination check report.
    """
    output:
        VAL_HTML = 'contamination_check/{input_format}/report/html-empty.iob'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/report',
        input_format = lambda wildcards: wildcards.input_format,
        input_type = config['input']['type']
    run:
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.quality_checks.htmlreportercontamination import HtmlReporterContamination

        # Suffix for hybrid data
        if params.input_type == 'hybrid':
            suffix_title = 'Illumina' if params.input_format == 'fastq_pe' else 'ONT'
        else:
            suffix_title = None
        section = HtmlReporterContamination.generate_empty_section(suffix_title)
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule contamination_check_dump_summary_info:
    """
    Dumps the summary information for the contamination check in tabular format.
    """
    input:
        INFORMS_kraken = rules.contamination_check_kraken2_run.output.INFORMS,
        INFORMS_species = rules.contamination_check_kraken_report_parser.output.INFORMS
    output:
        FILE = 'contamination_check/{input_format}/summary/summary_out.{ext}'
    params:
        dir_ = lambda wildcards: f'contamination_check/{wildcards.input_format}/summary',
        input_type = config['input']['type'],
        input_format = lambda wildcards: wildcards.input_format,
        allowed_species = config['contamination_check'].get('allowed_species', None),
        ext = lambda wildcards: wildcards.ext
    run:
        informs_kraken2 = snakemakeutils.load_object(Path(input.INFORMS_kraken))
        informs = snakemakeutils.load_object(Path(input.INFORMS_species))

        suffix = f'_{params.input_format}' if params.input_type == 'hybrid' else ''
        summary_data = [
            (f'kraken2{suffix}_expected_taxon', informs['expected'][0]),
            (f'kraken2{suffix}_expected_taxon_occurrence', informs['expected'][1]),
            (f'kraken2{suffix}_contaminants_warn', str(informs['contaminants_warn'])),
            (f'kraken2{suffix}_contaminants_fail', str(informs['contaminants_fail'])),
        ]
        if params.allowed_species is not None:
            summary_data.append((f'kraken2{suffix}_allowed', str(informs['allowed'])))
        summary_data.extend([
            (f'kraken2{suffix}_tool_version', informs_kraken2['_name_full']),
            (f'kraken2{suffix}_db', informs_kraken2['database']['name']),
            (f'kraken2{suffix}_last_update', informs_kraken2['database']['last_update'])
        ])
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), f'contamination_check{suffix}')
