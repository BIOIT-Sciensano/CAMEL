from pathlib import Path
from camel.app.snakemake import snakemakeutils
from camel.app.pipeline.step import Step


rule species_determination_analysis:
    """
    This rule is used to determine the species based on the cgST.
    """
    input:
        TSV_profile_matches = 'species_determination/input/tsv_profile_matches.io',
        TSV_taxonomic = 'species_determination/input/tsv_taxonomic.io'
    output:
        TSV = 'species_determination/tool/tsv.io',
        INFORMS = 'species_determination/tool/informs.io' # species_determination.OUTPUT_INFORMS
    params:
        dir_ = 'species_determination/tool'
    run:
        from camel.app.tools.pipelines.yersinia.speciesdetermination import SpeciesDetermination
        detector = SpeciesDetermination()
        detector.add_input_files({
            'TSV_profile_matches': snakemakeutils.load_object(Path(str(input.TSV_profile_matches))),
            'TSV_taxonomic': snakemakeutils.load_object(Path(str(input.TSV_taxonomic)))})
        step = Step(rule_name=str(rule), tool=detector, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(detector, output)

rule species_determination_report:
    """
    This rule is used to determine the species determination HTML report.
    """
    input:
        INFORMS_analysis = rules.species_determination_analysis.output.INFORMS,
        TSV_analysis = rules.species_determination_analysis.output.TSV
    output:
        VAL_HTML = 'species_determination/report/html.iob' # species_determination.OUTPUT_REPORT
    params:
        dir_ = 'species_determination/report'
    run:
        from camel.app.tools.pipelines.yersinia.speciesdeterminationreporter import SpeciesDeterminationReporter
        reporter = SpeciesDeterminationReporter()
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule species_determination_report_empty:
    """
    This rule is used to create an empty report for the species determination.
    """
    output:
        VAL_HTML = 'species_determination/report/html-empty.iob' # species_determination.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Species determination', Path(output.VAL_HTML))

rule species_determination_create_summary:
    """
    Creates the tabular summary output for the species determination.
    """
    input:
        INFORMS = rules.species_determination_analysis.output.INFORMS
    output:
        FILE = 'species_determination/summary/summary.{ext}' # species_determination.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))['best_match']
        data_summary = [(f'species_determination_{key}', str(v)) for key, v in informs.items()]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'species_determination')
