from pathlib import Path

import pandas as pd
import json

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly


rule resfinder4_run:
    """
    Runs the ResFinder4 tool.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA,
        DIR = config['resfinder4']['db']
    output:
        TSV_genes = 'resfinder4/tsv-genes.io',
        TSV_point = 'resfinder4/tsv-point.io',
        TSV_pheno_species = 'resfinder4/tsv-pheno-species.io',
        TSV_pheno_general = 'resfinder4/tsv-pheno-general.io',
        INFORMS = 'resfinder4/informs.io' # resfinder4.OUTPUT_INFORMS
    params:
        dir_ = 'resfinder4',
        species = config['resfinder4'].get('species'),
        point = config['resfinder4'].get('point', True),
        min_id = config['resfinder4'].get('min_identity'),
        min_cov = config['resfinder4'].get('min_cov')
    run:
        from camel.app.tools.resfinder.resfinder import ResFinder
        resfinder = ResFinder()
        snakemakeutils.add_pickle_input(resfinder, 'FASTA', Path(input.FASTA))
        resfinder.add_input_files({'DIR': [ToolIODirectory(Path(input.DIR))]})
        resfinder.update_parameters(min_cov=0.9, acquired=True, point=params.point)
        if params.species is not None:
            resfinder.update_parameters(species=f'"{params.species}"')
        if params.min_id is not None:
            resfinder.update_parameters(threshold=params.min_id)
        if params.min_cov is not None:
            resfinder.update_parameters(min_cov=params.min_cov)
        step = Step(rule_name=str(rule), tool=resfinder, dir_=Path(str(params.dir_)))
        step.run()
        if params.point:
            snakemakeutils.dump_tool_outputs(resfinder, output)
        else:
            snakemakeutils.dump_tool_outputs(resfinder, output,
                keys=[key for key in output.keys() if key not in ('TSV_point', 'TSV_pheno_species')])
            snakemakeutils.dump_object([], Path(output.TSV_point))
            snakemakeutils.dump_object([], Path(output.TSV_pheno_species))

rule resfinder4_reporter:
    """
    Creates the output report for ResFinder4.
    """
    input:
        TSV_genes = rules.resfinder4_run.output.TSV_genes,
        TSV_point = rules.resfinder4_run.output.TSV_point,
        TSV_pheno_species = rules.resfinder4_run.output.TSV_pheno_species,
        TSV_pheno_general = rules.resfinder4_run.output.TSV_pheno_general,
        INFORMS_resfinder = rules.resfinder4_run.output.INFORMS
    output:
        VAL_HTML = 'resfinder4/html.iob' # resfinder4.OUTPUT_REPORT
    params:
        dir_ = 'resfinder4',
        point= config['resfinder4'].get('point',True)
    run:
        from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
        reporter = ResFinderReporter()
        snakemakeutils.add_pickle_inputs(reporter, input, excluded_keys = None if params.point else ['TSV_point', 'TSV_pheno_species'])
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule resfinder4_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'resfinder4/html-empty.iob' # resfinder4.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('ResFinder4', Path(output.VAL_HTML))

rule resfinder4_create_summary:
    """
    Creates a tabular summary output for ResFinder4.
    """
    input:
        INFORMS = rules.resfinder4_run.output.INFORMS,
        TSV_genes = rules.resfinder4_run.output.TSV_genes,
        TSV_point = rules.resfinder4_run.output.TSV_point
    output:
        FILE = 'resfinder4/summary_resfinder.{ext}' # resfinder4.OUTPUT_SUMMARY
    params:
        point = config['resfinder4'].get('point', True),
        ext = lambda wildcards: wildcards.ext
    run:
        tsv_genes = snakemakeutils.load_object(Path(input.TSV_genes))[0].path
        if params.point:
            tsv_point = snakemakeutils.load_object(Path(input.TSV_point))[0].path
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_genes = pd.read_table(tsv_genes)

        # Format hits
        if params.ext == 'json':
            hits_data = data_genes.to_dict('records')
        elif params.ext == 'tsv':
            hits_data = json.dumps(data_genes.astype(str).values.tolist())
        else:
            raise ValueError(f'Invalid format: {params.ext}')

        # Summary output
        summary_data = [
            ('resfinder4_genes', ', '.join(list(data_genes['Resistance gene'])) if not data_genes.empty else '-'),
            ('resfinder4_genes_hits', hits_data),
            ('resfinder4_tool_version', informs['_name_full']),
            ('resfinder4_db_version_date', informs['db_version_resfinder']),
            ('resfinder4_db_version_name', informs['db_version_name'])
        ]
        if params.point:
            data_mutations = pd.read_table(tsv_point)
            if params.ext == 'tsv':
                summary_data.append(('resfinder4_mutations', ', '.join(list(data_mutations['Mutation'])) if not data_mutations.empty else '-'))
            elif params.ext == 'json':
                summary_data.append(('resfinder4_mutations', {mutation['Mutation']: mutation['Resistance'].split(', ') for _, mutation in data_mutations.iterrows()} if not data_mutations.empty else '-'))
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), 'resfinder4')
