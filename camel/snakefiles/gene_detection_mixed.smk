from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake.step import Step
from camel.snakefiles import gene_detection

rule gene_detection_mixed_aggregate:
    """
    TODO
    """
    input:
        HITS_kma = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method='kma'),
        HITS_blast = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method='blast'),
        INFORMS_db = gene_detection.OUTPUT_DB_INFORMS
    output:
        VAL_hits = 'gene_detection/{db}/mixed/hits.iob',  # gene_detection.OUTPUT_HITS_METHOD
        INFORMS = 'gene_detection/{db}/mixed/informs.io',
        VAL_agg_hits = 'gene_detection/{db}/mixed/agg-hits-summary.iob',
        TSV_agg_hits = 'gene_detection/{db}/mixed/tsv-agg-hits-summary.iob'
    params:
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/mixed',
        min_identity = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_coverage', 60),
    run:
        from camel.app.tools.pipelines.genedetection.mixedhitaggregating import MixedHitAggregating
        mixed_hit = MixedHitAggregating()
        snakemakeutils.add_pickle_inputs(mixed_hit, input)
        step = Step(rule_name=str(rule), tool=mixed_hit, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(mixed_hit, output)

# rule gene_detection_export_informs:
#     """
#     TODO
#     """
#     input:
#         INFORMS_kma = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method='kma'),
#         INFORMS_blast = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method='blast')
#     output:
#         VAL_hits = 'gene_detection/{db}/mixed/hits.iob',  # gene_detection.OUTPUT_HITS_METHOD
#         INFORMS = 'gene_detection/{db}/mixed/informs.io',
#         VAL_agg_hits = 'gene_detection/{db}/mixed/agg-hits-summary.iob',
#         TSV_agg_hits = 'gene_detection/{db}/mixed/tsv-agg-hits-summary.iob'
#     params:
#         dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/mixed',
#         min_identity = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_percent_identity', 90),
#         min_coverage = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_coverage', 60),
#     run:
#         from camel.app.tools.pipelines.genedetection.mixedhitaggregating import MixedHitAggregating
#         mixed_hit = MixedHitAggregating()
#         snakemakeutils.add_pickle_inputs(mixed_hit, input)
#         step = Step(rule_name=str(rule), tool=mixed_hit, dir_=Path(str(params.dir_)))
#         step.run()
#         snakemakeutils.dump_tool_outputs(mixed_hit, output)