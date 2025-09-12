from pathlib import Path

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly
from camel.scripts.salmonellapipeline.snakefile import serotyping_seqsero2


rule serotyping_seqsero2_run_fasta:
    """
    Runs SeqSero2 to perform serotyping on FASTA input.
    Note: The detection method is always (genome assembly) 'kmer' for FASTA input.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TXT = 'serotyping/seqsero2/kmer/txt.io',
        INFORMS = 'serotyping/seqsero2/kmer/informs.io' # serotyping_seqsero2.OUTPUT_KMER_INFORMS
    params:
        mode = 'kmer',
        dir_ = 'serotyping/seqsero2/kmer',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2()
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db_path_seqsero2))]})
        snakemakeutils.add_pickle_input(seqsero_tool, 'FASTA', Path(input.FASTA))
        seqsero_tool.update_parameters(mode=params.mode)
        step = Step(rule_name=str(rule), tool=seqsero_tool, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(seqsero_tool, output)

rule serotyping_seqsero2_run_fastq:
    """
    Runs SeqSero2 to perform serotyping on FASTQ input.
    Note: The detection method can be 'allele' or 'kmerread' for FASTQ input.
    """
    input:
        IO = 'fq_dict.io'
    output:
        TXT = 'serotyping/seqsero2/{mode}/txt.io',
        INFORMS = 'serotyping/seqsero2/{mode}/informs.io'
    params:
        mode = lambda wildcards: wildcards.mode,
        dir_ = lambda wildcards: f'serotyping/seqsero2/{wildcards.mode}',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2()
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db_path_seqsero2))]})
        if config['input_type'] == 'illumina':
            seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(Path(input.IO), key_pe='FASTQ_PE'))
        elif config['input_type'] == 'ont':
            seqsero_tool.add_input_files(SnakePipelineUtils.extracts_fq_input(
                Path(input.IO), key_se='FASTQ_ONT', read_type='SE'))
        seqsero_tool.update_parameters(mode=str(params.mode))
        step = Step(rule_name=str(rule), tool=seqsero_tool, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(seqsero_tool, output)

rule serotyping_seqsero2_dump_summary_info:
    """
    This rule creates a simple output summary for the SeqSero2 serotyping analysis.
    """
    input:
        TXT_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.TXT,
        INFORMS_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.INFORMS,
        TXT_seqsero2_allele = rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        INFORMS_seqsero2_allele = rules.serotyping_seqsero2_run_fastq.output.INFORMS.format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        TXT_seqsero2_kmerread = rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_seqsero2_kmerread = rules.serotyping_seqsero2_run_fastq.output.INFORMS.format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else []
    output:
        FILE = 'serotyping/seqsero2/summary/summary_seqsero2.{ext}' # serotyping_seqsero2.OUTPUT_SUMMARY
    params:
        dir_ = 'serotyping_seqsero2/summary',
        ext = lambda wildcards: wildcards.ext
    run:
        # Parse required SeqSero2 output
        informs_seqsero2_kmer = snakemakeutils.load_object(Path(input.INFORMS_seqsero2_kmer))
        tsv_results_full = serotyping_seqsero2.seqsero2_output_parser(snakemakeutils.load_object(Path(input.TXT_seqsero2_kmer))[0].path, 'seqsero2_kmer')

        # parse facultative SeqSero2 output
        if 'fasta' not in config['input']:
            files = []
        if input.TXT_seqsero2_allele:
            files = [(snakemakeutils.load_object(Path(str(input.TXT_seqsero2_allele)))[0].path, 'seqsero2_allele')]
            if input.TXT_seqsero2_kmerread:  # kmerread output can only be present if allele output is present
                files.append((snakemakeutils.load_object(Path(str(input.TXT_seqsero2_kmerread)))[0].path, 'seqsero2_kmerread'))
            for args_tuple in files:
                tsv_results = serotyping_seqsero2.seqsero2_output_parser(*args_tuple)
                tsv_results_full.extend(tsv_results)

        items = [tuple(item.split('\t')) for item in tsv_results_full]
        rows_out = [
            *items,
            ('seqsero2_tool_version', informs_seqsero2_kmer['_name']),
            ('seqsero2_db_version', informs_seqsero2_kmer['last_update_date'])
        ]
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'seqsero2')

rule serotyping_seqsero2_report:
    """
    Creates the HTML report for SeqSero2.
    If Illumina FASTQ data is available, all three modes are executed.
    If ONT FASTQ data is available, the allele mode is not executed because although the tool supposedly provides 
    support, it doesn't actually work.
    For FASTA and hybrid input, only the FASTA input is used.
    """
    input:
        TXT_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.TXT,
        TXT_seqsero2_allele = lambda wildcards: rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='allele') if config['input_type'] == 'illumina' else [],  # exclude fasta, ONT & hybrid
        TXT_seqsero2_kmerread = lambda wildcards: rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='kmerread') if config['input_type'] in ('ont', 'illumina') else [],
        INFORMS_serotyping_seqsero2 = rules.serotyping_seqsero2_run_fasta.output.INFORMS
    output:
        VAL_HTML = 'serotyping/seqsero2/report/html.iob' # serotyping_seqsero2.OUTPUT_REPORT
    params:
        dir_ = 'serotyping/seqsero2/report',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path']
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter

        reporter = SeqSero2Reporter()
        reporter.add_input_files({'DIR_seqsero2': [ToolIODirectory(Path(params.db_path_seqsero2))]})
        snakemakeutils.add_pickle_inputs(reporter, input, excluded_keys=['TXT_seqsero2_allele', 'TXT_seqsero2_kmerread'])
        if input.TXT_seqsero2_allele:
            snakemakeutils.add_pickle_input(reporter, 'TXT_seqsero2_allele', Path(str(input.TXT_seqsero2_allele)))
        if input.TXT_seqsero2_kmerread:
            snakemakeutils.add_pickle_input(reporter, 'TXT_seqsero2_kmerread', Path(str(input.TXT_seqsero2_kmerread)))
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule serotyping_seqsero2_report_empty:
    """
    Creates an empty HTML report for the SeqSero2 serotyping analysis.
    """
    output:
        VAL_HTML = 'serotyping/seqsero2/report/html-empty.iob' # serotyping_seqsero2.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter
        SnakePipelineUtils.create_empty_report_section(SeqSero2Reporter.TITLE, Path(output.VAL_HTML))
