from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.viralconsensuspipeline.snakefile import dbcomparison, iterativemapping, refselection


checkpoint db_comparison_split_consensus:
    """
    Splits the final consensus sequence by segment.
    """
    input:
        FASTA = iterativemapping.OUTPUT_FASTA_CONSENSUS_FINAL_TRIMMED
    output:
        DIR = directory('db_comparison/consensus')
    params:
        segments_db = refselection.get_segments(Path(config['ref_selection']['db'])) if config['ref_selection'].get(
            'db') is not None else []
    run:
        from Bio import SeqIO
        from camel.app.loggers import logger

        # Save each segment present in the database to a separate file
        Path(output.DIR).mkdir(parents=True)
        path_fasta = snakemakeutils.load_object(Path(input.FASTA))[0].path
        with path_fasta.open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                segment = seq.id.split('-')[-1]
                if segment not in params.segments_db:
                    logger.warning(f'Consensus segment not found in the database: {segment}')
                    continue
                if set(str(seq.seq).upper()) <= {'N'}:
                    logger.warning(f'Consensus segment without called bases: {segment}')
                    continue
                SeqIO.write(seq, Path(output.DIR) / f'{segment}.fasta', 'fasta')

rule db_comparison_mash_dist:
    """
    Calculates the mash distance between the consensus segment and all database sequences of that segment.
    """
    input:
        FASTA = 'db_comparison/consensus/{segment}.fasta',
        DB = lambda wildcards: dbcomparison.get_db_path(config, wildcards.segment, 'mash')
    output:
        TSV = 'db_comparison/mash_dist/{segment}/tsv.io',
        INFORMS = 'db_comparison/mash_dist/{segment}/informs.io' # dbcomparison.OUTPUT_INFORMS_MASH
    params:
        dir_ = lambda wildcards: f'db_comparison/mash_dist/{wildcards.segment}'
    threads: 4
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.mash.mashdist import MashDist
        mash_dist = MashDist()
        mash_dist.add_input_files({'FASTA': [ToolIOFile(Path(input.FASTA))], 'DB': [ToolIOFile(Path(str(input.DB)))]})
        mash_dist.update_parameters(threads=threads)
        step = Step(rule_name=str(rule), tool=mash_dist, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(mash_dist, output)

rule db_comparison_select_candidates:
    """
    Selects the database sequences with the lowest mash distance as candidates for the alignment.
    """
    input:
        TSV = rules.db_comparison_mash_dist.output.TSV,
        INFORMS = rules.db_comparison_mash_dist.output.INFORMS
    output:
        TXT = 'db_comparison/candidates/{segment}/candidates.txt'
    params:
        max_candidates = 100
    run:
        import pandas as pd
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        cols = snakemakeutils.load_object(Path(input.INFORMS))['cols']
        data_mash = pd.read_table(path_tsv, names=cols)
        data_mash.nsmallest(params.max_candidates, 'distance')['reference_id'].to_csv(
            output.TXT, index=False, header=False)

rule db_comparison_extract_candidates:
    """
    Extracts the candidate sequences from the database.
    """
    input:
        TXT = rules.db_comparison_select_candidates.output.TXT,
        FASTA = lambda wildcards: dbcomparison.get_db_path(config, wildcards.segment, 'fasta')
    output:
        FASTA = 'db_comparison/candidates/{segment}/fasta.io'
    params:
        dir_ = lambda wildcards: f'db_comparison/candidates/{wildcards.segment}'
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
        samtools_faidx = SamtoolsFastaIndex()
        samtools_faidx.add_input_files({
            'FASTA': [ToolIOFile(Path(str(input.FASTA)))],
            'TXT': [ToolIOFile(Path(input.TXT))]
        })
        samtools_faidx.update_parameters(output_filename='candidates.fasta', symlink_input=False)
        step = Step(rule_name=str(rule), tool=samtools_faidx, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(samtools_faidx, output)

rule db_comparison_align:
    """
    Aligns the candidate sequences to the consensus segment using minimap2.
    """
    input:
        FASTA_REF = 'db_comparison/consensus/{segment}.fasta',
        FASTA = rules.db_comparison_extract_candidates.output.FASTA
    output:
        PAF = 'db_comparison/align/{segment}/paf.io',
        INFORMS = 'db_comparison/align/{segment}/informs.io' # dbcomparison.OUTPUT_INFORMS_ALIGN
    params:
        dir_ = lambda wildcards: f'db_comparison/align/{wildcards.segment}'
    threads: 4
    run:
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.minimap2.minimap2alignment import Minimap2Alignment
        minimap2 = Minimap2Alignment()
        minimap2.add_input_files({'FASTA_REF': [ToolIOFile(Path(input.FASTA_REF))]})
        snakemakeutils.add_io_inputs(minimap2, input, keys=['FASTA'])
        minimap2.update_parameters(cs_tag=True, threads=threads)
        step = Step(rule_name=str(rule), tool=minimap2, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(minimap2, output)

rule db_comparison_compare:
    """
    Ranks the database sequences by the number of differences with the consensus sequence.
    """
    input:
        PAF = lambda wildcards: [
            str(rules.db_comparison_align.output.PAF).format(segment=segment)
            for segment in dbcomparison.get_segments(checkpoints)],
        FASTA = lambda wildcards: [
            f'db_comparison/consensus/{segment_}.fasta' for segment_ in dbcomparison.get_segments(checkpoints)],
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        TSV = 'db_comparison/compare/tsv.io',
        TSV_equal = 'db_comparison/compare/tsv-equal.io',
        INFORMS = 'db_comparison/compare/informs.io'
    params:
        dir_ = 'db_comparison/compare'
    run:
        import itertools
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.pipelines.viral_consensus.dbcomparison import DbComparison
        db_comparison = DbComparison()
        db_comparison.add_input_files({
            'PAF': list(itertools.chain(*[snakemakeutils.load_object(Path(p)) for p in input.PAF])),
            'FASTA': [ToolIOFile(Path(p)) for p in input.FASTA],
            'DB': [ToolIODirectory(Path(input.DB))]
        })
        step = Step(rule_name=str(rule), tool=db_comparison, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(db_comparison, output)

rule db_comparison_report:
    """
    Creates the report section for the database comparison.
    """
    input:
        TSV = rules.db_comparison_compare.output.TSV,
        TSV_equal = rules.db_comparison_compare.output.TSV_equal,
        INFORMS_db_comparison = rules.db_comparison_compare.output.INFORMS,
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        VAL_HTML = 'db_comparison/report/html.iob' # dbcomparison.OUTPUT_REPORT
    params:
        dir_ = 'db_comparison/report'
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.viral_consensus.reporterdbcomparison import ReporterDbComparison
        reporter = ReporterDbComparison()
        snakemakeutils.add_io_inputs(reporter, input, excluded_keys=['DB'])
        reporter.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule db_comparison_report_empty:
    """
    Creates an empty report section when the database comparison is disabled.
    """
    output:
        VAL_HTML = 'db_comparison/report/html-empty.iob' # dbcomparison.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Database comparison', Path(output.VAL_HTML))

rule db_comparison_summary:
    """
    Creates the summary output for the database comparison.
    """
    input:
        TSV = rules.db_comparison_compare.output.TSV,
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        FILE = 'db_comparison/summary/summary.{ext}' # dbcomparison.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import pandas as pd
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        data_matches = pd.read_table(path_tsv, keep_default_na=False, na_values=[''])  # Segment 'NA' is not missing
        data_summary = [('db_comparison_database', Path(input.DB).name)]
        for segment, data_seg in data_matches.groupby('segment'):
            closest = data_seg.iloc[0]
            data_summary.extend([
                (f'db_comparison-{segment}', closest['accession']),
                (f'db_comparison_identity-{segment}', round(float(closest['identity']), 2)),
                (f'db_comparison_differences-{segment}', int(closest['differences'])),
                (f'db_comparison_nb_equal-{segment}', int((data_seg['rank'] == 1).sum()))
            ])
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'db_comparison')
