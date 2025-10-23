from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule confindr_run:
    """
    Runs ConFindr on the input FASTQ datasets.
    """
    input:
        IO_FASTQ = 'fq_dict.io'
    output:
        CSV = 'confindr/tool/csv.io',
        INFORMS = 'confindr/tool/informs.io' # confindr.OUTPUT_INFORMS
    params:
        dir_ = 'confindr/tool',
        db = config.get('confindr', {}).get('db'),
        input_type = config['input_type']
    threads: 4
    run:
        from camel.app.scriptutils.fastqinput import FastqInput
        from camel.app.tools.confindr.confindr import ConFindr

        confindr_ = ConFindr()
        confindr_.update_parameters(rmlst=True, threads=threads)

        # Add database
        if params.db is None:
            raise RuntimeError('ConFindr database not specified in config')
        path_db = Path(params.db)
        if not path_db.exists():
            raise FileNotFoundError(f'ConFindr database not found: {path_db}')
        confindr_.update_parameters(rmlst=True, databases=str(path_db))

        # Add input files
        if params.input_type in ('hybrid', 'illumina'):
            fq_in = FastqInput.from_fq_dict(Path(input.IO_FASTQ), 'illumina')
            confindr_.add_input_files({'FASTQ_PE': fq_in.pe})
        else:
            fq_in = FastqInput.from_fq_dict(Path(input.IO_FASTQ),'ont')
            confindr_.add_input_files({'FASTQ_SE': fq_in.se})
            confindr_.update_parameters(data_type='Nanopore',quality_cutoff=12)

        # Run the tool
        step = Step(rule_name=str(rule), tool=confindr_, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(confindr_, output)

rule confindr_report:
    """
    Creates an output report for the ConFindr analysis.
    """
    input:
        INFORMS_confindr = rules.confindr_run.output.INFORMS
    output:
        HTML = 'confindr/report/html.iob' # confindr.OUTPUT_REPORT
    params:
        dir_ = 'confindr/report'
    run:
        from camel.app.tools.confindr.confindrreporter import ConFindrReporter
        reporter = ConFindrReporter()
        reporter.update_parameters(input_type=config['input_type'])
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        snakemakeutils.add_pickle_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule confindr_report_empty:
    """
    Creates an empty ConFindr report when the analysis is disabled.
    """
    output:
        HTML = 'confindr/report/html-empty.iob' # confindr.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('ConFindr', Path(output.HTML), 2)

rule confindr_create_summary:
    """
    Creates the tabular summary output for the ConFindr assay.
    """
    input:
        INFORMS_confindr = rules.confindr_run.output.INFORMS
    output:
        FILE = 'confindr/summary/summary_confindr.{ext}' # confindr.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        # Collect informs
        informs = snakemakeutils.load_object(Path(input.INFORMS_confindr))
        rows_out = [
            ('genus', informs['Genus']),
            ('nb_contam_snvs', informs['NumContamSNVs']),
            ('contam_status', informs['ContamStatus']),
            ('tool_version', informs['_name_full']),
            ('db_version', informs['DatabaseDownloadDate'])
        ]
        rows_out = [(f'confindr_{k}', v) for k, v in rows_out]
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'confindr')
