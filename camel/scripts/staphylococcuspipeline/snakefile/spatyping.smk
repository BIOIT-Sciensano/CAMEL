from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import assembly

rule spa_typing_blastn:
    """
    Runs BLASTN to align the sequences against the spa typing database.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA,
        DB_BLAST = config['spa_typing']['db']
    output:
        TSV = 'spa_typing/blastn/tsv.io',
        INFORMS = 'spa_typing/blastn/informs.io' # 'spatyping_workflow.OUTPUT_INFORMS
    params:
        dir_ = 'spa_typing/blastn'
    run:
        from camel.app.tools.blast.blastn import Blastn
        from camel.app.tools.spatyping.spatyping import SpaTyping
        from camel.app.core.io.tooliofile import ToolIOFile
        blastn = Blastn()
        snakemakeutils.add_pickle_input(blastn, 'FASTA', Path(input.FASTA))
        blastn.add_input_files({'DB_BLAST': [ToolIOFile(Path(input.DB_BLAST))]})
        blastn.update_parameters(
            output_format=SpaTyping.BLASTN_OUTPUT_FORMAT,
            num_alignments=100_000,
            dust='no',
            task='blastn')
        step = Step(rule_name=str(rule), tool=blastn, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(blastn, output)

rule spa_typing_run:
    """
    Determines the spa-type based on the blastn output.
    """
    input:
        TSV = rules.spa_typing_blastn.output.TSV,
        CSV_profiles = config['spa_typing']['profiles']
    output:
        VAL_hits = 'spa_typing/detection/hits.iob',
        INFORMS = 'spa_typing/detection/informs.io'
    params:
        dir_ = 'spa_typing/detection'
    run:
        from camel.app.core.io.tooliofile import ToolIOFile
        from camel.app.tools.spatyping.spatyping import SpaTyping
        spatyping = SpaTyping()
        snakemakeutils.add_pickle_input(spatyping, 'TSV', Path(input.TSV))
        spatyping.add_input_files({'CSV_profiles': [ToolIOFile(Path(input.CSV_profiles))]})
        step = Step(rule_name=str(rule), tool=spatyping, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(spatyping, output)

rule spa_typing_report:
    """
    Creates a report for the spa-typing assay.
    """
    input:
        VAL_hits = rules.spa_typing_run.output.VAL_hits,
        INFORMS_spa_typing = rules.spa_typing_run.output.INFORMS
    output:
        VAL_HTML = 'spa_typing/report/html.iob' # spatyping_workflow.OUTPUT_REPORT
    params:
        dir_ = 'spa_typing/report'
    run:
        from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter
        reporter = SpaTypingReporter()
        snakemakeutils.add_pickle_inputs(reporter, input)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule spa_typing_report_empty:
    """
    Creates an empty report when spatyping is disabled.
    """
    output:
        VAL_HTML = 'spa_typing/report/html-empty.iob' # spatyping_workflow.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('<i>spa</i> typing', Path(output.VAL_HTML))

rule spa_typing_summary:
    """
    Exports the summary information for the spa-typing assay
    """
    input:
        INFORMS = rules.spa_typing_run.output.INFORMS
    output:
        FILE = 'spa_typing/summary_spatyping.{ext}' # spatyping_workflow.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS))
        data_summary = [
            ('spa_type', informs['spa_type'] if informs['spa_type'] is not None else '-'),
            ('spa_type_repeats', informs['spa_type_repeats'] if informs['spa_type'] is not None else '-'),
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'spatyping')
