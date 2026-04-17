from pathlib import Path

from pandas.errors import EmptyDataError

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import mobsuite


rule mobsuite_mob_recon:
    """
    Runs the MOB-recon tool.
    """
    input:
        FASTA = mobsuite.INPUT_FASTA, # mobsuite.INPUT_FASTA
        DB = config['mob_suite']['db']
    output:
        TSV = 'mob_suite/tool/tsv.io',
        TSV_contigs = 'mob_suite/tool/tsv-contigs.io',
        FASTA = 'mob_suite/tool/fasta.io',
        INFORMS = 'mob_suite/tool/informs.io' # mobsuite.OUTPUT_INFORMS
    params:
        dir_ = 'mob_suite/tool'
    threads: 4
    run:
        from camel.app.tools.mobsuite.mobrecon import MOBRecon
        mob_recon = MOBRecon()
        snakemakeutils.add_io_input(mob_recon,'FASTA', Path(input.FASTA))
        mob_recon.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        mob_recon.update_parameters(num_threads=threads)
        step = Step(rule_name=str(rule), tool=mob_recon, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(mob_recon, output)

rule mobsuite_mob_recon_reporter:
    """
    Creates the output report for MOB-recon.
    """
    input:
        TSV = rules.mobsuite_mob_recon.output.TSV,
        TSV_contigs = rules.mobsuite_mob_recon.output.TSV_contigs,
        FASTA = rules.mobsuite_mob_recon.output.FASTA,
        INFORMS_mob_recon = rules.mobsuite_mob_recon.output.INFORMS
    output:
        HTML = 'mob_suite/report/html.iob' # mobsuite.OUTPUT_REPORT
    params:
        dir_ = 'mob_suite/report',
        contig_report = config.get('mob_suite', {}).get('contig_report', False)
    run:
        from camel.app.tools.mobsuite.mobreconreporter import MOBReconReporter
        reporter = MOBReconReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        reporter.update_parameters(contig_report=params.contig_report)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule mobsuite_mob_recon_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'mob_suite/report/html-empty.iob' # mobsuite.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('MOB-recon', Path(output.VAL_HTML))

rule mobsuite_create_summary:
    """
    Creates a tabular summary output for MOB-suite.
    """
    input:
        TSV = rules.mobsuite_mob_recon.output.TSV,
        INFORMS = rules.mobsuite_mob_recon.output.INFORMS
    output:
        FILE = 'mob_suite/summary/summary_mob_suite.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import re
        import pandas as pd
        from camel.app.tools.mobsuite.mobreconreporter import MOBReconReporter

        # Parse TSV output
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        try:
            data_mobsuite = pd.read_table(path_tsv)
            if 'primary_cluster_id' not in data_mobsuite.columns:
                raise EmptyDataError
            primary_cluster_ids = list(data_mobsuite['primary_cluster_id'])
        except EmptyDataError:
            logger.info(f'No plasmids detected by MOB-suite')
            primary_cluster_ids = []

        # Parse informs
        informs_in = snakemakeutils.load_object(Path(input.INFORMS))

        # Create summary output
        rows_out = [
            ('mob_suite_primary_cluster_ids', ', '.join(primary_cluster_ids) if len(primary_cluster_ids) > 0 else '-' ),
            ('mob_suite_predicted_plasmid_contigs', ', '.join(ctg for ctg, status in informs_in['contig_report'].items() if status is not None)),
            ('mob_suite_tool_version', informs_in['_name_full'])
        ]
        if primary_cluster_ids and params.ext == 'json':
            data_mobsuite['id'] = data_mobsuite['sample_id'].apply(lambda x: re.search('.*:(.*)',x).group(1))
            data_mobsuite['id'] = data_mobsuite['id'].apply(lambda x: MOBReconReporter.format_plasmid_id(x))
            cluster_info = [{col:d.get('fmt', lambda x: x)(row[col]) for col, d in MOBReconReporter.COLUMN_MAPPING.items()} for row in data_mobsuite.to_dict('records')]
            rows_out.append(('mob_suite_overview', cluster_info))
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'mob_suite')

rule mobsuite_report_genomic_context:
    """
    Reports the genomic context for detected genes.
    """
    input:
        TSV = 'mob_suite/genomic_context/input/tsv.io',
        INFORMS = 'mob_suite/genomic_context/input/informs.io',
        INFORMS_mob_recon = rules.mobsuite_mob_recon.output.INFORMS
    output:
        HTML = 'mob_suite/genomic_context/html.iob' # mobsuite.OUTPUT_CONTEXT_REPORT
    params:
        dir_ = 'mob_suite/genomic_context',
        detection_method = config['gene_detection']['options']['method'],
        input_type = config.get('input_type', 'illumina')
    run:
        from camel.app.tools.mobsuite.genomiccontext import GenomicContext

        genomic_context = GenomicContext()
        genomic_context.add_input_files({'TSV': snakemakeutils.load_object(Path(input.TSV))})
        genomic_context.add_input_informs({
            'mob_recon': snakemakeutils.load_object(Path(input.INFORMS_mob_recon)),
            'dbs': snakemakeutils.load_object(Path(input.INFORMS))
        })
        genomic_context.update_parameters(
            detection_method=str(params.detection_method), input_type=str(params.input_type))
        step = Step(rule_name=str(rule), tool=genomic_context, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(genomic_context, output)

rule mobsuite_report_genomic_context_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'mob_suite/genomic_context/html-empty.iob' # mobsuite.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Genomic context', Path(output.VAL_HTML), 2)
