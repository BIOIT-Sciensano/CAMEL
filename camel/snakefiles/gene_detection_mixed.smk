from pathlib import Path

from camel.app.toolkits.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection

# Include workflows for the different detection methods
# include: gene_detection.SNAKEFILE
include: gene_detection.SNAKEFILE_KMA
include: gene_detection.SNAKEFILE_BLAST

rule gene_detection_db_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA = 'gene_detection/{db}/db_manager/fasta.io',
        FASTA_clustered = 'gene_detection/{db}/db_manager/fasta-clust.io',
        INFORMS = 'gene_detection/{db}/db_manager/informs.iob' # gene_detection.OUTPUT_DB_INFORMS
    params:
        db_path = lambda wildcards: config['gene_detection']['dbs'][wildcards.db]['path'],
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}'
    run:
        from camel.app.core.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        db_manager = DBManager()
        db_manager.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
        step = Step(
            rule_name=str(rule),
            tool=db_manager,
            dir_=Path(str(params.dir_)),
            wildcards=wildcards)
        step.run()
        snakemakeutils.dump_tool_outputs(db_manager, output)

rule gene_detection_get_hits_mixed:
    """
    Retrieves the hits based on the detection method species in the configuration.
    """
    input:
        INFORMS_DB = rules.gene_detection_db_manager.output.INFORMS,
        VAL_hits = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db)),
        INFORMS_hits = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db))
    output:
        VAL_hits = 'gene_detection/{db}/{method}/hit_selection/hits-standardized.iob',
        INFORMS = 'gene_detection/{db}/{method}/hit_selection/informs.io' # gene_detection.OUTPUT_INFORMS
    run:
        import shutil
        shutil.copyfile(str(input.VAL_hits), output.VAL_hits)
        # Add a tag for the database to distinguish commands in the output
        informs = snakemakeutils.load_object(Path(str(input.INFORMS_hits)))
        informs['_tag'] = snakemakeutils.load_object(Path(str(input.INFORMS_DB)))['title']
        snakemakeutils.dump_object(informs, Path(output.INFORMS))

rule gene_detection_map_names_mixed:
    """
    Maps the standardized names (seq_X) to the original ones.
    Renames the TSV file based on the sample name and database.
    """
    input:
        INFORMS_db = rules.gene_detection_db_manager.output.INFORMS,
        HITS = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db),
    output:
        HITS = 'gene_detection/{db}/{method}/hit_selection/selected-hits.iob', # gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS,
        TSV = 'gene_detection/{db}/{method}/metadata/tsv.io'
    params:
        dir_working = lambda wildcards: f'gene_detection/{wildcards.db}/{wildcards.method}/metadata',
        sample_name = config['input']['sample_name'],
        db_config = lambda wildcards: config['gene_detection']['dbs'][wildcards.db]
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.core.io.tooliovalue import ToolIOValue
        from camel.app.core.utils import fileutils
        informs_db = snakemakeutils.load_object(Path(input.INFORMS_db))
        hits = snakemakeutils.load_object(Path(input.HITS))

        # Map standardized names to original ones
        hits_updated = []
        for hit in [io.value for io in hits]:
            seq_id = hit.locus
            hit.locus = informs_db['mapping'].get_metadata(seq_id, 'allele')
            hit.accession = informs_db['mapping'].get_metadata(seq_id, 'accession', '-')
            if params.db_config.get('metadata') is not None:
                key = params.db_config['metadata']['key']
                hit.add_metadata(params.db_config['metadata']['name'], informs_db['mapping'].get_metadata(seq_id, key))
            hits_updated.append(hit)
        snakemakeutils.dump_object([ToolIOValue(hit) for hit in hits_updated], Path(output.HITS))

        # Save tabular output
        if len(hits_updated) >= 1:
            output_path = Path(str(params.dir_working)) / 'hits-{}-{}.tsv'.format(
                fileutils.make_valid(params.sample_name), fileutils.make_valid(informs_db['name']))
            with output_path.open('w') as handle:
                handle.write('\t'.join(hits_updated[0].table_column_names))
                handle.write('\n')
                for hit in hits_updated:
                    handle.write('\t'.join(hit.to_table_row()))
                    handle.write('\n')
            snakemakeutils.dump_object([ToolIOFile(output_path)], Path(output.TSV))
        else:
            snakemakeutils.dump_object([], Path(output.TSV))

rule gene_detection_get_column_names_mixed:
    """
    Retrieves the column names for the gene detection output.
    This method is necessary in case output needs to be generated when no hits are detected.
    """
    output:
        INFORMS_columns = 'gene_detection/{db}/{method}/report/informs-columns.io' # gene_detection.OUTPUT_COLUMNS
    params:
        detection_method = lambda wildcards: wildcards.method,
        db_config = lambda wildcards: config['gene_detection']['dbs'][wildcards.db]
    run:
        from camel.app.toolkits.blast.blasthitstatistics import BlastHitStatistics
        from camel.app.toolkits.genedetection.genedetectionblasthit import GeneDetectionBlastHit
        from camel.app.toolkits.genedetection.genedetectionkmahit import GeneDetectionKMAHit

        # Create empty hit
        if params.detection_method == 'blast':
            empty_hit = GeneDetectionBlastHit('Locus', None, BlastHitStatistics('subject', 0, '', 'query', 0, 0, '', 0.0, '+'))
        elif params.detection_method == 'kma':
            empty_hit = GeneDetectionKMAHit('DB_cluster', 'Locus', None, '', 0, 0, 0.0, 0.0, 0.0)
        else:
            raise ValueError(f"Invalid detection method: {params.detection_method}")

        # Add metadata columns
        if 'metadata' in params.db_config:
            empty_hit.add_metadata(params.db_config['metadata']['name'], '')

        # Save column names
        columns = empty_hit.html_column_names
        snakemakeutils.dump_object(columns, Path(output.INFORMS_columns))

rule reporter_mixed:
    input:
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS,
        INFORMS_columns_k = 'gene_detection/{db}/report/kma/informs-columns.io',
        INFORMS_columns_b = 'gene_detection/{db}/report/blast/informs-columns.io',
        HITS_k = 'gene_detection/{db}/hit_selection/kma/selected-hits.iob',
        HITS_b = 'gene_detection/{db}/hit_selection/blast/selected-hits.iob',
        TSV_k='gene_detection/{db}/metadata/kma/tsv.io',
        TSV_b='gene_detection/{db}/metadata/blast/tsv.io'
    output:
        VAL_HTML='gene_detection/{db}/report/mixed_report_html.iob'
    params:
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/report',
        config_data = lambda wildcards: config['gene_detection']['dbs'][wildcards.db],
        input_type = config['input']['type'],
        detection_method = lambda wildcards: GeneDetectionUtils.get_detection_method_key(config, wildcards.db)
    run:
        from camel.app.core.reports.htmlreportsection import HtmlReportSection
        from camel.app.core.io.tooliovalue import ToolIOValue

        db_info = snakemakeutils.load_object(Path(input.INFORMS_db_info))
        section = HtmlReportSection(db_info['title'],3)
        section.add_paragraph('Analysis for mixed genes is included in the report.')
        snakemakeutils.dump_object([ToolIOValue(section)],Path(output.VAL_HTML))

        # from camel.app.tools.pipelines.genedetection.htmlreportergenedetection import HtmlReporterGeneDetection
        #
        # # Create step
        # reporter = HtmlReporterGeneDetection()
        # step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)), wildcards=wildcards)
        #
        # # Parameters
        # if 'force_detection_method' in params.config_data:
        #     reporter.update_parameters(forced_detection_method = params.config_data['force_detection_method'])
        # if params.config_data.get('hidden', False) is True:
        #     reporter.update_parameters(hidden=True)
        # if params.input_type == 'fasta' and params.detection_method in ('kma',):
        #     reporter.update_parameters(pseudo_reads=True)
        #
        # # Optional message
        # if 'message' in params.config_data:
        #     reporter.update_parameters(message=params.config_data['message']['content'])
        # if ('message' in params.config_data) and ('category' in params.config_data['message']):
        #     reporter.update_parameters(message_category=params.config_data['message']['category'])
        #
        # # Run tool
        # snakemakeutils.add_pickle_inputs(reporter, input)
        # step.run()
        # snakemakeutils.dump_tool_outputs(reporter, output)
