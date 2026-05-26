from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.snakefiles import sequence_typing, assembly
from camel.app.core.snakemake import snakemakeutils


rule serogroup_determination_analysis:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        hits = expand(sequence_typing.OUTPUT_HITS,
            locus_type='DNA',
            detection_method = config['sequence_typing']['options']['method'],
            scheme = sorted([k for k in config['sequence_typing']['dbs'].keys() if k.startswith('serogroup')]))
    output:
        INFORMS = 'serogroup_determination/legacy/tool/informs.io'
    params:
        serogroups = sorted([k for k in config['sequence_typing']['dbs'].keys() if k.startswith('serogroup')]),
        dir_ = 'serogroup_determination/legacy/tool'
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdetermination import SerogroupDetermination
        detector = SerogroupDetermination()
        for hits_output, serogroup in zip(input.hits, params.serogroups):
            detector.add_input_files({f'hits_{serogroup}': snakemakeutils.load_object(Path(hits_output))})
        step = Step(rule_name=str(rule), tool=detector, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(detector, output)

rule serogroup_determination_report:
    """
    This rule is used to determine the serogroup based on the detected genes.
    """
    input:
        INFORMS_analysis = rules.serogroup_determination_analysis.output.INFORMS
    output:
        VAL_HTML = 'serogroup_determination/legacy/report/html.iob' # serogroup_determination.OUTPUT_LEGACY_REPORT
    run:
        from camel.app.tools.pipelines.neisseria.serogroupdeterminationreporter import SerogroupDeterminationReporter
        reporter = SerogroupDeterminationReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule serogroup_determination_report_empty:
    """
    This rule is used to create an empty report for the serogroup determination
    """
    output:
        VAL_HTML = 'serogroup_determination/legacy/report/html-empty.iob' # serogroup_determination.OUTPUT_LEGACY_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Serogroup determination (legacy)', Path(output.VAL_HTML))

rule serogroup_capsule_tool:
    """
    Runs the serogroup capsule tool.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TSV = 'serogroup_determination/capsule/tool/tsv.io',
        JSON = 'serogroup_determination/capsule/tool/json.io',
        INFORMS = 'serogroup_determination/capsule/tool/informs.io' # serogroup_determination.OUTPUT_INFORMS
    run:
        from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsule import CharacterizeNeisseriaCapsule
        capsule_typer = CharacterizeNeisseriaCapsule()
        snakemakeutils.add_io_inputs(capsule_typer, input)
        step = Step(rule_name=str(rule), tool=capsule_typer, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(capsule_typer, output)

rule serogroup_capsule_tool_report:
    """
    Creates an output report for the capsule tool.
    """
    input:
        TSV = rules.serogroup_capsule_tool.output.TSV,
        JSON = rules.serogroup_capsule_tool.output.JSON,
        INFORMS_detector = rules.serogroup_capsule_tool.output.INFORMS
    output:
        HTML = 'serogroup_determination/capsule/report/html.iob' # serogroup_determination.OUTPUT_REPORT
    run:
        from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsulereporter import \
            CharacterizeNeisseriaCapsuleReporter
        reporter = CharacterizeNeisseriaCapsuleReporter()
        snakemakeutils.add_io_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule serogroup_capsule_tool_report_empty:
    """
    Creates an empty report when the analysis is disabled.
    """
    output:
        VAL_HTML = 'serogroup_determination/capsule/report/html-empty.iob' # serogroup_determination.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Capsule characterization', Path(output.VAL_HTML))

rule serogroup_determination_summary:
    """
    Collects the summary information for the serogroup determination.
    """
    input:
        INFORMS_capsule = rules.serogroup_capsule_tool.output.INFORMS,
        INFORMS_legacy = rules.serogroup_determination_analysis.output.INFORMS
    output:
        FILE = 'serogroup_determination/summary/summary_out.{ext}' # serogroup_determination.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        rows_out = []

        # Capsule tool
        informs_capsule = snakemakeutils.load_object(Path(input.INFORMS_capsule))
        rows_out.extend([
            ('serogroup_capsule', informs_capsule['detected_serogroup']),
            ('serogroup_capsule_genes', informs_capsule['genes_present'])
        ])

        # Legacy detection
        informs_legacy = snakemakeutils.load_object(Path(input.INFORMS_legacy))
        rows_out.extend([
            ('serogroup_legacy', informs_legacy['detected_serogroup']),
            ('serogroup_legacy_nb_hits', informs_legacy['serogroups_sorted'][0]['nb_hits']),
            ('serogroup_legacy_nb_hits_perfect', informs_legacy['serogroups_sorted'][0]['nb_hits_perfect']),
            ('serogroup_legacy_total_loci', informs_legacy['serogroups_sorted'][0]['nb_loci_total'])
        ])

        snakemakeutils.export_summary(rows_out, Path(output.FILE), str(params.ext), 'serogroup')
