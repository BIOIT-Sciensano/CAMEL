from pathlib import Path

import pandas as pd
import json

from camel.app.core.io.tooliodirectory import ToolIODirectory
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import amrfinder


rule amrfinder_run:
    """
    Runs the AMRFinder tool.
    """
    input:
        FASTA = amrfinder.INPUT_FASTA,
        DIR = config['amrfinder']['db']
    output:
        TSV = 'amrfinder/tool/tsv.io', # amrfinder.OUTPUT_TSV
        INFORMS = 'amrfinder/tool/informs.io' # amrfinder.OUTPUT_INFORMS
    params:
        dir_ = 'amrfinder/tool',
        organism = config['amrfinder'].get('species')
    run:
        from camel.app.tools.amrfinder.amrfinder import AMRFinder
        amrfinder = AMRFinder()
        snakemakeutils.add_io_input(amrfinder,'FASTA', Path(input.FASTA))
        amrfinder.add_input_files({'DIR': [ToolIODirectory(Path(input.DIR))]})
        if params.organism is not None:
            amrfinder.update_parameters(organism=params.organism)
        step = Step(rule_name=str(rule), tool=amrfinder, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(amrfinder, output)

rule amrfinder_reporter:
    """
    Creates the output report for AMRFinder.
    """
    input:
        TSV = rules.amrfinder_run.output.TSV,
        INFORMS_amrfinder = rules.amrfinder_run.output.INFORMS
    output:
        HTML = 'amrfinder/report/html.iob' # amrfinder.OUTPUT_REPORT
    params:
        dir_ = 'amrfinder/report'
    run:
        from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter
        reporter = AMRFinderReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule amrfinder_report_empty:
    """
    Creates an empty report when this analysis is disabled.
    """
    output:
        VAL_HTML = 'amrfinder/html-empty.iob' # amrfinder.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('AMRFinder', Path(output.VAL_HTML))

rule amrfinder_dump_summary_info:
    """
    Dumps the summary information for the ResFinder workflow in tabular format.
    """
    input:
        TSV = rules.amrfinder_run.output.TSV,
        INFORMS = rules.amrfinder_run.output.INFORMS
    output:
        FILE = 'amrfinder/summary_amrfinder.{ext}' # amrfinder.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        data_amr = pd.read_table(path_tsv)

        # Extract the informs
        informs = snakemakeutils.load_object(Path(input.INFORMS))

        # Parse perfect & other hits
        if not data_amr.empty:
            data_amr['is_perfect'] = data_amr.apply(
                lambda x: x['% Coverage of reference'] == 100.0 and x['% Identity to reference'] == 100.0, axis=1)
            hits_perfect = list(data_amr[data_amr['is_perfect']]['Element symbol']) if sum(data_amr['is_perfect']) > 0 else []
            hits_other = list(data_amr[~data_amr['is_perfect']]['Element symbol']) if sum(~data_amr['is_perfect']) > 0 else []
        else:
            hits_perfect = []
            hits_other = []

        # Format hits
        if params.ext == 'tsv':
            data_hits = json.dumps(data_amr.drop(columns=['is_perfect']).astype(str).values.tolist()) if not data_amr.empty else '-'
        elif params.ext == 'json':
            data_hits = list(data_amr.to_dict('records')) if not data_amr.empty else []
        else:
            raise ValueError(f'Invalid ext: {params.ext}')

        # Summary output
        data = [
            ('amr_hits_perfect', ', '.join(hits_perfect) if len(hits_perfect) > 0 else '-'),
            ('amr_hits_other', ', '.join(hits_other) if len(hits_other) > 0 else '-'),
            ('amr_genes_hits', data_hits),
            ('amr_tool_version', informs['_name_full']),
            ('amr_db_version', informs['db_version'])
        ]
        snakemakeutils.export_summary(data, Path(output.FILE), str(params.ext), 'amrfinder')
