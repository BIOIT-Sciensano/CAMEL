from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.shigellapipeline.snakefile import shigatyper


rule shigatyper_run:
    """
    Runs the ShigaTyper assay.
    """
    input:
        IO = shigatyper.get_input(config)
    output:
        TSV = 'shigatyper/tool/tsv.io',
        TSV_HITS = 'shigatyper/tool/tsv_hits.io',
        INFORMS = 'shigatyper/tool/informs.io' # shigatyper.OUTPUT_INFORMS
    params:
        dir_ = 'shigatyper/tool',
        input_type = config['input']['type']
    run:
        from camel.app.tools.pipelines.shigella.shigatyper import ShigaTyper
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput

        typer = ShigaTyper()

        # Extract FASTQ paths to add them as ShigaTyper input
        if params.input_type == 'illumina':
            fq_in = FastqInput.from_fq_dict(Path(input.IO), 'illumina')
            typer.add_input_files({'FASTQ_PE': fq_in.pe})
        elif params.input_type == 'ont':
            fq_in = FastqInput.from_fq_dict(Path(input.IO),'ont')
            typer.update_parameters(ont=True)
            typer.add_input_files({'FASTQ_SE': fq_in.se})
        else:
            # When the input type is FASTA the simulated reads are used as input
            snakemakeutils.add_io_input(typer,'FASTQ_PE', Path(input.IO))

        # Run tool
        step = Step(rule_name=str(rule), tool=typer, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(typer, output)

rule shigatyper_report:
    """
    Creates an output report for the ShigaTyper analysis.
    """
    input:
        TSV = rules.shigatyper_run.output.TSV,
        TSV_HITS = rules.shigatyper_run.output.TSV_HITS,
        INFORMS_shigatyper = rules.shigatyper_run.output.INFORMS
    output:
        HTML = 'shigatyper/report/html.iob' # shigatyper.OUTPUT_REPORT
    params:
        dir_ = 'shigatyper/report',
        input_type = config['input']['type']
    run:
        from camel.app.tools.pipelines.shigella.shigatyperreporter import ShigaTyperReporter
        reporter = ShigaTyperReporter()
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        snakemakeutils.add_io_inputs(reporter, input)
        if params.input_type == 'fasta':
            reporter.update_parameters(pseudo_reads=True)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule shigatyper_report_empty:
    """
    Creates an empty ShigaTyper report when the analysis is disabled.
    """
    output:
        HTML = 'shigatyper/report/html-empty.iob' # shigatyper.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('ShigaTyper', Path(output.HTML), 3)

rule shigatyper_create_summary:
    """
    Creates the tabular summary output for the ShigaTyper assay.
    """
    input:
        INFORMS_shigatyper = rules.shigatyper_run.output.INFORMS
    output:
        FILE = 'shigatyper/summary_shigatyper.{ext}' # shigatyper.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        informs = snakemakeutils.load_object(Path(input.INFORMS_shigatyper))
        data_summary = [
            ('shigatyper_prediction', str(informs['species'])),
            ('shigatyper_tool_version', informs['_name_full']),
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'shigatyper')
