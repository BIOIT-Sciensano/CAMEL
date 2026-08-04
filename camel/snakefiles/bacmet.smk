from pathlib import Path

import pandas as pd

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly

rule bacmet_pickle_db:
    """
    Creates a pickle object for the database.
    """
    input:
        DB = config['bacmet']['db']
    output:
        DB = 'bacmet/db.io'
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        snakemakeutils.dump_object([ToolIODirectory(Path(input.DB))], Path(output.DB))

rule bacmet_prodigal:
    """
    Runs Prodigal to predict CDS.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        FASTA = 'bacmet/prodigal/tool/fasta.io',
        INFORMS = 'bacmet/prodigal/tool/informs.io' # bacmet.OUTPUT_PRODIGAL_INFORMS
    run:
        from camel.app.tools.prodigal.prodigal import Prodigal
        prodigal = Prodigal()
        snakemakeutils.add_io_inputs(prodigal, input)
        step = Step(rule_name=str(rule), tool=prodigal, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(prodigal, output)

rule bacmet_prodigal_report:
    """
    Creates an output report for Prodigal.
    """
    input:
        FASTA = rules.bacmet_prodigal.output.FASTA,
        INFORMS_prodigal = rules.bacmet_prodigal.output.INFORMS
    output:
        HTML = 'bacmet/prodigal/report/html.iob' # bacmet.OUTPUT_PRODIGAL_REPORT
    run:
        from camel.app.tools.prodigal.prodigalreporter import ProdigalReporter
        reporter = ProdigalReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule bacmet_prodigal_report_empty:
    """
    Creates an empty output report for Prodigal.
    """
    output:
        HTML = 'bacmet/prodigal/report/html-empty.iob' # bacmet.OUTPUT_PRODIGAL_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Prodigal', Path(output.HTML))

rule bacmet_blastp:
    """
    Runs blastp to identify protein matches.
    """
    input:
        FASTA = rules.bacmet_prodigal.output.FASTA,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        TSV = 'bacmet/blastp/tsv.io',
        INFORMS = 'bacmet/blastp/informs.io' # bacmet.OUTPUT_INFORMS
    params:
        fmt = '6 pident sseqid sseq slen qseqid qstart qend'
    threads: 4
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.blast.blastp import Blastp

        # Create & run tool
        blastp = Blastp()
        blastp.update_parameters(output_format=f'"{params.fmt}"', threads=threads)
        snakemakeutils.add_io_input(blastp,'FASTA', Path(input.FASTA))
        path_db = next(snakemakeutils.load_object(Path(input.DB))[0].path.glob('*.fasta'))
        blastp.add_input_files({'DB_BLAST': [ToolIOFile(path_db)]})
        step = Step(rule_name=str(rule), tool=blastp, dir_=snakemakeutils.get_rule_dir(output))
        step.run()

        # Dump output
        snakemakeutils.dump_io_outputs(blastp, output)

rule bacmet_filter_blastp:
    """
    Filters the hits detected by blastp.
    """
    input:
        TSV = rules.bacmet_blastp.output.TSV,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        TSV = 'bacmet/hit_filtering/tsv.io',
        INFORMS = 'bacmet/hit_filtering/informs.io'
    params:
        dir_ = 'bacmet/hit_filtering',
        min_cov = 90,
        min_id = 75,
        cols = 'pident sseqid sseq slen qseqid qstart qend'
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        tsv_in = snakemakeutils.load_object(Path(input.TSV))[0].path
        data_in = pd.read_table(tsv_in, names=params.cols.split(' '))
        data_in['perc_covered'] = data_in.apply(lambda x: 100.0 * float(len(x['sseq'])) / x['slen'], axis=1)
        data_in_filt = data_in[data_in['pident'] > params.min_id].copy()
        data_in_filt = data_in_filt[data_in_filt['perc_covered'] > params.min_cov].copy()
        data_in_filt['BacMet_ID'] = data_in_filt['sseqid'].apply(lambda x: x.split('|')[0])

        # Parse metadata
        dir_db = snakemakeutils.load_object(Path(input.DB))[0].path
        tsv_meta = next(dir_db.glob('*.tsv'))
        data_meta = pd.read_table(tsv_meta)
        data_out = pd.merge(data_in_filt, data_meta, on='BacMet_ID')
        path_out = Path(params.dir_, 'bacmet.tsv')
        data_out.to_csv(path_out, sep='\t', index=False)
        snakemakeutils.dump_object([ToolIOFile(path_out)], Path(output.TSV))
        snakemakeutils.dump_object(
            {'params': {'min_cov': params.min_cov, 'min_id': params.min_id}}, Path(output.INFORMS))

rule bacmet_report:
    """
    Creates the output report for BacMet. 
    """
    input:
        TSV = rules.bacmet_filter_blastp.output.TSV,
        INFORMS_blastp = rules.bacmet_blastp.output.INFORMS,
        INFORMS_filtering = rules.bacmet_filter_blastp.output.INFORMS,
        DB = rules.bacmet_pickle_db.output.DB
    output:
        HTML = 'bacmet/report/html.iob' # bacmet.OUTPUT_REPORT
    run:
        from camel.app.tools.pipelines.klebsiella.bacmetreporter import BacMetReporter
        reporter = BacMetReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule bacmet_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'bacmet/report/html-empty.iob' # bacmet.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('BacMet', Path(output.VAL_HTML))

rule bacmet_create_summary:
    """
    Creates a tabular summary output file for the BacMet assay.
    """
    input:
        TSV = rules.bacmet_filter_blastp.output.TSV
    output:
        TSV = 'bacmet/summary_bacmet.tsv' # bacmet.OUTPUT_SUMMARY
    run:
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        data_in = pd.read_table(path_tsv)
        with open(output.TSV, 'w') as handle:
            handle.write('\t'.join(['bacmet_genes', ', '.join(sorted(list(data_in['Gene_name'])))]))
            handle.write('\n')
