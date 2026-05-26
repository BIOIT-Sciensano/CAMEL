from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import gene_detection


rule gene_detection_blast_blastn:
    """
    Aligns sequences agains the database with blastn.
    """
    input:
        FASTA = gene_detection.INPUT_FASTA,
        DB_BLAST = 'gene_detection/{db}/db_manager/fasta-clust.io'
    output:
        ASN = 'gene_detection/{db}/blastn/asn.io',
        INFORMS = 'gene_detection/{db}/blastn/informs.io'
    params:
        task = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('blastn', {}).get('task', 'megablast'),
        blast_reads = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('blastn', {}).get('blast_reads', False)
    run:
        from camel.app.tools.blast.blastn import Blastn
        blastn = Blastn()
        snakemakeutils.add_io_inputs(blastn, input)
        step = Step(rule_name=str(rule), tool=blastn, dir_=snakemakeutils.get_rule_dir(output))
        blastn.update_parameters(threads=1, task=str(params.task), max_target_seqs=1 if params.blast_reads else 20_000)
        step.run()
        blastn.informs['Task'] = params.task
        snakemakeutils.dump_io_outputs(blastn, output)

rule gene_detection_blast_tsv_generation:
    """
    Generates tabular output format to extract hit statistics.
    """
    input:
        ASN = rules.gene_detection_blast_blastn.output.ASN
    output:
        TSV = 'gene_detection/{db}/tsv_generation/tsv.io'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter()
        snakemakeutils.add_io_inputs(blast_formatter, input)
        step = Step(rule_name=str(rule), tool=blast_formatter, dir_=snakemakeutils.get_rule_dir(output))
        blast_formatter.update_parameters(output_format='"7 pident sseqid sseq slen qseqid qstart qend score"')
        step.run()
        snakemakeutils.dump_io_outputs(blast_formatter, output)

rule gene_detection_blast_hit_filtering:
    """
    Filters hits based on percent identity and query coverage.
    Extracts the hit information based on the database metadata.
    """
    input:
        TSV = rules.gene_detection_blast_tsv_generation.output.TSV,
        INFORMS_blastn = rules.gene_detection_blast_blastn.output.INFORMS
    output:
        VAL_Hits = 'gene_detection/{db}/hit_filtering/blast-hits.iob',
        INFORMS = 'gene_detection/{db}/blast/informs.io' # gene_detection.OUTPUT_INFORMS_METHOD
    params:
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/hit_filtering',
        min_percent_identity = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('blastn', {}).get('min_percent_identity', 90),
        min_coverage = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('blastn', {}).get('min_coverage', 60),
        filtering_method = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params',{}).get('blastn',{}).get('filtering_method', 'cluster'),
        score_nb_of_hits = lambda wildcards: config['gene_detection']['dbs'][wildcards.db].get('params', {}).get('blastn', {}).get('score_nb_of_hits', 5)
    run:
        from camel.app.tools.pipelines.genedetection.blasthitfiltering import BlastHitFiltering
        hit_filtering = BlastHitFiltering()
        snakemakeutils.add_io_inputs(hit_filtering, input)
        step = Step(rule_name=str(rule), tool=hit_filtering, dir_=Path(str(params.dir_)))

        # Update parameters
        hit_filtering.update_parameters(
            min_percent_identity=str(params.min_percent_identity),
            min_coverage=str(params.min_coverage),
            filtering_method=str(params.filtering_method),
            score_nb_of_hits=str(params.score_nb_of_hits)
        )
        # Run tool
        step.run()
        snakemakeutils.dump_io_outputs(hit_filtering, output)

        # Add the informs from the filtering to the existing ones with the blastn command
        informs = snakemakeutils.load_object(Path(input.INFORMS_blastn))
        for key, value in hit_filtering.informs.items():
            if key.startswith('_'):
                continue
            informs[key] = value
        snakemakeutils.dump_object(informs, Path(output.INFORMS))

rule gene_detection_blast_text_alignment_generation:
    """
    Generates alignments in the text format.
    """
    input:
        ASN = rules.gene_detection_blast_blastn.output.ASN
    output:
        TXT = 'gene_detection/{db}/alignment_generation/txt.io'
    run:
        from camel.app.tools.blast.blastformatter import BlastFormatter
        blast_formatter = BlastFormatter()
        snakemakeutils.add_io_inputs(blast_formatter, input)
        step = Step(rule_name=str(rule), tool=blast_formatter, dir_=snakemakeutils.get_rule_dir(output))
        blast_formatter.update_parameters(output_format='0', num_alignments=20000)
        step.run()
        snakemakeutils.dump_io_outputs(blast_formatter, output)

rule gene_detection_blast_text_alignment_extraction:
    """
    Extracts a text alignment for the selected hits and attaches them to the hit objects.
    """
    input:
        TXT = rules.gene_detection_blast_text_alignment_generation.output.TXT,
        VAL_Hits = rules.gene_detection_blast_hit_filtering.output.VAL_Hits
    output:
        VAL_Hits =  'gene_detection/{db}/blast/hits.iob' # gene_detection.OUTPUT_HITS_METHOD
    params:
        dir_ = lambda wildcards: f'gene_detection/{wildcards.db}/alignment_generation'
    run:
        from camel.app.tools.pipelines.genedetection.alignmentextractor import AlignmentExtractor
        alignment_extractor = AlignmentExtractor()
        snakemakeutils.add_io_inputs(alignment_extractor, input)
        step = Step(rule_name=str(rule), tool=alignment_extractor, dir_=Path(str(params.dir_)))
        step.run()
        hits_with_alignment = []
        for io_value, alignment in zip(
                snakemakeutils.load_object(Path(input.VAL_Hits)),
                alignment_extractor.tool_outputs['TXT']
        ):
            io_value.value.alignment_path = Path(alignment.path)
            hits_with_alignment.append(io_value)
        snakemakeutils.dump_object(hits_with_alignment, Path(output.VAL_Hits))
