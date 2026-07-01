from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake.step import Step
from camel.snakefiles import gene_detection

rule gene_detection_mixed_aggregate:
    """
    Aggregate the KMA and Blast hits into a single mixed hit object.
    """
    input:
        HITS_kma = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method='kma'),
        HITS_blast = lambda wildcards: gene_detection.OUTPUT_HITS_METHOD.format(db=wildcards.db, method='blast'),
        INFORMS_db = gene_detection.OUTPUT_DB_INFORMS
    output:
        VAL_hits = 'gene_detection/{db}/mixed/hits.iob',  # gene_detection.OUTPUT_HITS_METHOD
        INFORMS = 'gene_detection/{db}/mixed/informs_pre.io',
        VAL_hits_complete = 'gene_detection/{db}/mixed/hits-complete.iob'  # gene_detection.OUTPUT_HITS_COMPLETE
    params:
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/mixed',
        min_identity = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('mixed', {}).get('min_coverage', 60)
    run:
        from camel.app.tools.pipelines.genedetection.mixedhitaggregating import MixedHitAggregating
        mixed_hit = MixedHitAggregating()
        snakemakeutils.add_pickle_inputs(mixed_hit, input)
        step = Step(rule_name=str(rule), tool=mixed_hit, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(mixed_hit, output)

rule gene_detection_mixed_export_informs:
    """
    Adds the commands informs so that they are displayed correctly in the report.
    """
    input:
        INFORMS = rules.gene_detection_mixed_aggregate.output.INFORMS,
        INFORMS_kma = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method='kma'),
        INFORMS_blast = lambda wildcards: gene_detection.OUTPUT_INFORMS_METHOD.format(db=wildcards.db, method='blast')
    output:
        INFORMS = 'gene_detection/{db}/mixed/informs.io'  # gene_detection.OUTPUT_INFORMS_METHOD
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        informs_kma = snakemakeutils.load_object(Path(str(input.INFORMS_kma)))
        informs_blast = snakemakeutils.load_object(Path(str(input.INFORMS_blast)))
        informs['_kma'] = snakemakeutils.load_object(Path(str(input.INFORMS_kma)))
        informs['_blast'] = snakemakeutils.load_object(Path(str(input.INFORMS_blast)))
        informs['_command'] = '\n'.join(i['_command'] for i in (informs_kma, informs_blast))
        snakemakeutils.dump_object(informs, Path(output.INFORMS))
