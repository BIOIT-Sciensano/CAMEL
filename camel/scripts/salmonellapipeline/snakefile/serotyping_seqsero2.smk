from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.snakefiles import assembly, read_simulation


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
        snakemakeutils.add_io_input(seqsero_tool,'FASTA', Path(input.FASTA))
        seqsero_tool.update_parameters(mode=params.mode)
        step = Step(rule_name=str(rule), tool=seqsero_tool, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(seqsero_tool, output)

def get_fastq_input(wildcards, input_type: str) -> tuple[str, bool]:
    """
    Returns the FASTQ input for SeqSero2.
    :param wildcards: Snakemake wildcards
    :param input_type: Input type
    :return: Path to input
    """
    if input_type == 'ont':
        if wildcards.mode == 'kmerread':
            return 'fq_dict.io', False
        return read_simulation.OUTPUT_FASTQ, True

    if input_type == 'fasta':
        return read_simulation.OUTPUT_FASTQ, True

    # Default case (e.g., illumina)
    return 'fq_dict.io', False

rule serotyping_seqsero2_run_fastq:
    """
    Runs SeqSero2 to perform serotyping on FASTQ input.
    Note: The detection method can be 'allele' or 'kmerread' for FASTQ input.
    """
    input:
        IO = lambda wildcards: get_fastq_input(wildcards, config['input']['type'])[0]
    output:
        TXT = 'serotyping/seqsero2/{mode}/txt.io',
        INFORMS = 'serotyping/seqsero2/{mode}/informs.io'
    params:
        mode = lambda wildcards: wildcards.mode,
        dir_ = lambda wildcards: f'serotyping/seqsero2/{wildcards.mode}',
        db_path_seqsero2 = config['serotyping']['seqsero2']['path'],
        input_type = config['input']['type'],
        reads_simulated = lambda wildcards: get_fastq_input(wildcards, config['input']['type'])[1]
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2

        seqsero_tool = SeqSero2()
        seqsero_tool.add_input_files({'DIR': [ToolIODirectory(Path(params.db_path_seqsero2))]})

        # Add FASTQ input
        input_path = Path(str(input.IO))
        if params.input_type == 'illumina':
            seqsero_tool.add_input_files(snakepipelineutils.extract_fq_input(input_path, key_pe='FASTQ_PE'))
        elif params.reads_simulated:
            seqsero_tool.add_input_files({'FASTQ_PE': snakemakeutils.load_object(input_path)})
        else:
            seqsero_tool.add_input_files(snakepipelineutils.extract_fq_input(
                input_path, key_se='FASTQ_ONT', read_type='SE'))

        # Update parameters and run tool
        seqsero_tool.update_parameters(mode=str(params.mode))
        step = Step(rule_name=str(rule), tool=seqsero_tool, dir_=Path(str(params.dir_)))
        step.run()

        # Update and store the informs
        seqsero_tool.informs['simulated'] = params.reads_simulated
        snakemakeutils.dump_io_outputs(seqsero_tool, output)

rule serotyping_seqsero2_dump_summary_info:
    """
    Creates the summary output for SeqSero2.
    """
    input:
        TXT_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.TXT,
        INFORMS_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.INFORMS,
        TXT_seqsero2_allele = rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='allele'),
        TXT_seqsero2_kmerread = rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='kmerread') if config['input']['type'] in ('ont', 'illumina') else []
    output:
        FILE = 'serotyping/seqsero2/summary/summary_seqsero2.{ext}' # serotyping_seqsero2.OUTPUT_SUMMARY
    params:
        dir_ = 'serotyping_seqsero2/summary',
        ext = lambda wildcards: wildcards.ext
    run:
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter

        # Parse required SeqSero2 output
        summary_data = []
        for key_file, path in input.items():
            if not key_file.startswith('TXT') or len(path) == 0:
                continue
            path_tsv = snakemakeutils.load_object(Path(path))[0].path
            seqsero2_out = SeqSero2Reporter.parse_seqsero_output(path_tsv)
            mode = key_file.split('_')[-1]
            for column in SeqSero2Reporter.COLS:
                col_key = column['key']
                short_name = column['name_short'].lower()
                summary_data.append((f"seqsero2_{mode}_{short_name}", seqsero2_out.get(col_key, '')))
            summary_data.append((f'seqsero2_{mode}_note', seqsero2_out.get('Note', 'n/a')))

        # Add tool information
        informs = snakemakeutils.load_object(Path(input.INFORMS_seqsero2_kmer))
        summary_data.extend([
            ('seqsero2_tool_version', informs['_name_full']),
            ('seqsero2_db_version', informs['last_update_date'])
        ])
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), 'seqsero2')

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
        INFORMS_seqsero2_kmer = rules.serotyping_seqsero2_run_fasta.output.INFORMS,
        TXT_seqsero2_allele = lambda wildcards: rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='allele'),
        INFORMS_seqsero2_allele = lambda wildcards: rules.serotyping_seqsero2_run_fastq.output.INFORMS.format(mode='allele'),
        TXT_seqsero2_kmerread = lambda wildcards: rules.serotyping_seqsero2_run_fastq.output.TXT.format(mode='kmerread') if config['input']['type'] in ('ont', 'illumina') else [],
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
        snakemakeutils.add_io_inputs(reporter, input, optionals=['TXT_seqsero2_kmerread', 'INFORMS_seqsero2_kmerread'])
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule serotyping_seqsero2_report_empty:
    """
    Creates an empty HTML report for the SeqSero2 serotyping analysis.
    """
    output:
        VAL_HTML = 'serotyping/seqsero2/report/html-empty.iob' # serotyping_seqsero2.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.pipelines.salmonella.seqsero2reporter import SeqSero2Reporter
        snakepipelineutils.create_empty_report_section(SeqSero2Reporter.TITLE, Path(output.VAL_HTML))
