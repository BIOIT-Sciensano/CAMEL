import itertools
import json
from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.viralconsensuspipeline.snakefile import nextclade3


rule antivirals_check_mutations:
    """
    Queries the detected antiviral mutations against the database.
    """
    input:
        TSV = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, 'TSV', config),
        INFORMS_subtype = lambda wildcards: nextclade3.get_informs_subtype(wildcards, checkpoints)
    output:
        JSON = 'antivirals/check/json.io'
    params:
        dir_ = 'antivirals/check',
        species = config.get('antivirals', {}).get('species'),
        db = config.get('antivirals', {}).get('db')
    run:
        from camel.app.tools.pipelines.viral_consensus.antiviralsdetection import AntiviralsDetection
        detection = AntiviralsDetection()
        detection.add_input_files({
            'DB': [ToolIODirectory(Path(params.db))],
            'TSV': list(itertools.chain(*[snakemakeutils.load_object(Path(tsv)) for tsv in input.TSV]))
        })

        # Check the subtypes
        informs_subtype = snakemakeutils.load_object(Path(str(input.INFORMS_subtype)))
        subtype = informs_subtype.get('subtype')
        if subtype is None:
            raise RuntimeError(f'A subtype should be detected to detected antiviral mutations.')
        detection.update_parameters(species=params.species, subtype=subtype)
        step = Step(rule_name=str(rule), tool=detection, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(detection, output)

rule antivirals_report:
    """
    Creates a HTML report for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        VAL_HTML = 'antivirals/report/html.iob' # antivirals.OUTPUT_REPORT
    params:
        dir_ = 'antivirals/report'
    run:
        from camel.app.tools.pipelines.viral_consensus.antiviralsreporter import AntiviralsReporter
        reporter = AntiviralsReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule antivirals_report_empty:
    """
    Creates an empty output report when the assay is disabled.
    """
    output:
        VAL_HTML = 'antivirals/report/html-empty.iob' # antivirals.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Antiviral resistance', Path(output.VAL_HTML))

rule antivirals_summary:
    """
    Creates the summary output for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        FILE = 'antivirals/summary/summary.{ext}' # antivirals.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        path_json = snakemakeutils.load_object(Path(input.JSON))[0].path
        with open(path_json) as handle:
            data_mutations = json.load(handle)
        data_summary = [
            ('antivirals_mutations', json.dumps(data_mutations['mutations'])),
            ('antivirals_associations', json.dumps(data_mutations['associations'])),
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'antivirals')
