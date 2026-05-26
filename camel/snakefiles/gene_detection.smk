from pathlib import Path

from camel.app.toolkits.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection

# Include workflows for the different detection methods
include: gene_detection.SNAKEFILE_BLAST
include: gene_detection.SNAKEFILE_KMA


# Common rules
rule gene_detection_db_manager:
    """
    Retrieves the FASTA file and the metadata from a database folder.
    """
    output:
        FASTA = 'gene_detection/{db}/db_manager/fasta.io',
        FASTA_clustered = 'gene_detection/{db}/db_manager/fasta-clust.io',
        INFORMS = 'gene_detection/{db}/db_manager/informs.iob' # gene_detection.OUTPUT_DB_INFORMS
    params:
        db_name = lambda wildcards: wildcards.db,
        db_path = lambda wildcards: config['gene_detection']['dbs'][wildcards.db]['path'],
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}'
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager

        if Path(str(params.db_path)).exists():
            db_manager = DBManager()
            db_manager.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db_path)))]})
            step = Step(
                rule_name=str(rule),
                tool=db_manager,
                dir_=Path(str(params.dir_)),
                wildcards=wildcards)
            step.run()
            snakemakeutils.dump_io_outputs(db_manager,output)
        else:
            snakemakeutils.dump_object([], Path(output.FASTA))
            snakemakeutils.dump_object([], Path(output.FASTA_clustered))
            snakemakeutils.dump_object({
                'name': params.db_name,
                'title': params.db_name,
            }, Path(output.INFORMS))


rule gene_detection_get_hits:
    """
    Retrieves the hits based on the detection method species in the configuration.
    """
    input:
        INFORMS_DB = rules.gene_detection_db_manager.output.INFORMS,
        VAL_hits = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db)),
        INFORMS_hits = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method=GeneDetectionUtils.get_detection_method_key(config, wildcards.db))
    output:
        VAL_hits = 'gene_detection/{db}/hit_selection/hits-standardized.iob',
        INFORMS = 'gene_detection/{db}/hit_selection/informs.io' # gene_detection.OUTPUT_INFORMS
    run:
        import shutil
        shutil.copyfile(str(input.VAL_hits), output.VAL_hits)
        # Add a tag for the database to distinguish commands in the output
        informs = snakemakeutils.load_object(Path(str(input.INFORMS_hits)))
        informs['_tag'] = snakemakeutils.load_object(Path(str(input.INFORMS_DB)))['title']
        snakemakeutils.dump_object(informs, Path(output.INFORMS))

rule gene_detection_map_names:
    """
    Maps the standardized names (seq_X) to the original ones.
    Renames the TSV file based on the sample name and database.
    """
    input:
        INFORMS_db = rules.gene_detection_db_manager.output.INFORMS,
        HITS = rules.gene_detection_get_hits.output.VAL_hits
    output:
        HITS = 'gene_detection/{db}/hit_selection/selected-hits.iob', # gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS,
        TSV = 'gene_detection/{db}/metadata/tsv.io'
    params:
        dir_working = lambda wildcards: f'gene_detection/{wildcards.db}/metadata',
        sample_name = config['input']['sample_name'],
        db_config = lambda wildcards: config['gene_detection']['dbs'][wildcards.db]
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camelcore.app.io.tooliovalue import ToolIOValue
        from camelcore.app.utils import fileutils
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

rule gene_detection_get_column_names:
    """
    Retrieves the column names for the gene detection output.
    This method is necessary in case output needs to be generated when no hits are detected.
    """
    output:
        INFORMS_columns = 'gene_detection/{db}/report/informs-columns.io' # gene_detection.OUTPUT_COLUMNS
    params:
        detection_method = lambda wildcards: GeneDetectionUtils.get_detection_method_key(config, wildcards.db),
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


rule gene_detection_report:
    """
    Creates HTML reports for the gene detection.
    """
    input:
        VAL_Hits = rules.gene_detection_map_names.output.HITS,
        TSV = rules.gene_detection_map_names.output.TSV,
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS,
        INFORMS_detection = rules.gene_detection_get_hits.output.INFORMS
    output:
        VAL_HTML = 'gene_detection/{db}/report/html.iob' # gene_detection.OUTPUT_REPORT
    params:
        config_data = lambda wildcards: config['gene_detection']['dbs'][wildcards.db],
        input_type = config['input']['type'],
        detection_method = lambda wildcards: GeneDetectionUtils.get_detection_method_key(config, wildcards.db)
    run:
        from camel.app.tools.pipelines.genedetection.htmlreportergenedetection import HtmlReporterGeneDetection

        # Create step
        reporter = HtmlReporterGeneDetection()
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output), wildcards=wildcards)

        # Parameters
        if 'force_detection_method' in params.config_data:
            reporter.update_parameters(forced_detection_method = params.config_data['force_detection_method'])
        if params.config_data.get('hidden', False) is True:
            reporter.update_parameters(hidden=True)
        if params.input_type == 'fasta' and params.detection_method in ('kma',):
            reporter.update_parameters(pseudo_reads=True)

        # Optional message
        if 'message' in params.config_data:
            reporter.update_parameters(message=params.config_data['message']['content'])
        if ('message' in params.config_data) and ('category' in params.config_data['message']):
            reporter.update_parameters(message_category=params.config_data['message']['category'])

        # Run tool
        snakemakeutils.add_io_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule gene_detection_create_empty_report:
    """
    Creates an empty HTML report for the gene detection.
    """
    input:
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS
    output:
        VAL_HTML = 'gene_detection/{db}/report/html-empty.iob' # gene_detection.OUTPUT_REPORT_EMPTY
    run:
        from camelcore.app.reports.htmlreportsection import HtmlReportSection
        from camelcore.app.io.tooliovalue import ToolIOValue
        db_info = snakemakeutils.load_object(Path(input.INFORMS_db_info))
        section = HtmlReportSection(db_info['title'], 3)
        section.add_paragraph('Analysis disabled')
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output.VAL_HTML))

rule gene_detection_dump_summary_info:
    """
    Dumps the summary information from the gene detection in tabular format.
    """
    input:
        INFORMS_hits = 'gene_detection/{db}/hit_selection/selected-hits.iob',
        INFORMS_db_info = rules.gene_detection_db_manager.output.INFORMS
    output:
        FILE = 'gene_detection/{db}/report/summary_out.{ext}' # gene_detection.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext,
        db_key = lambda wildcards: wildcards.db
    run:
        import json

        # Parse input informs
        informs_hits = snakemakeutils.load_object(Path(input.INFORMS_hits))
        informs_db_info = snakemakeutils.load_object(Path(input.INFORMS_db_info))

        # Collect detected hits
        if params.ext == 'tsv':
            # Tabular output -> List with metric values as a JSON string
            data_hits = json.dumps([hit.value.to_table_row() for hit in informs_hits])
        elif params.ext == 'json':
            # JSON output -> List of dictionaries with hit statistics
            data_hits = [hit.value.to_dict() for hit in informs_hits]
        else:
            raise ValueError(f'Invalid extension: {params.ext}')

        # Dump the output
        prepend_key = f'{params.db_key}_' if params.ext == 'tsv' else ''
        rows_out = [
            (f'{prepend_key}loci', data_hits),
            (f'{prepend_key}db_version', informs_db_info.get('last_change', informs_db_info.get('last_updated', 'n/a')))
        ]
        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), str(params.db_key))
