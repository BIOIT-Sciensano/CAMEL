from pathlib import Path

from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliovalue import ToolIOValue
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import mykrobe


rule mykrobe_run:
    """
    Runs the Mykrobe assay.
    """
    input:
        IO = mykrobe.get_input(config)
    output:
        CSV = 'mykrobe/tool/csv.io',
        INFORMS = 'mykrobe/tool/informs.io' # mykrobe.OUTPUT_INFORMS
    params:
        species = Path(config['mykrobe']['species']),
        db_dir = Path(config['mykrobe']['db']),
        input_type = config['input']['type']
    run:
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput
        from camel.app.tools.mykrobe.mykrobe import Mykrobe

        typer = Mykrobe()

        # Extract FASTQ paths to add them as Mykrobe input
        if params.input_type == 'illumina':
            fq_in = FastqInput.from_fq_dict(Path(input.IO), params.input_type)
            typer.add_input_files({'FASTQ_PE': fq_in.pe})
        if params.input_type == 'ont':
            fq_in = FastqInput.from_fq_dict(Path(input.IO), params.input_type)
            typer.add_input_files({'FASTQ_SE': fq_in.se})
        if params.input_type == 'fasta':
            snakemakeutils.add_io_input(typer,'FASTA', Path(input.IO))
        typer.add_input_files({
            'DIR': [ToolIODirectory(params.db_dir)],
            'SPECIES': [ToolIOValue(params.species)]
        })

        # Run tool
        step = Step(rule_name=str(rule), tool=typer, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_output(typer,'CSV', Path(output.CSV))
        # Informs can contain numpy objects -> sanitize first
        snakemakeutils.dump_object(snakemakeutils.sanitize_numpy(typer.informs), Path(output.INFORMS))

rule mykrobe_report:
    """
    Creates an output report for the Mykrobe analysis.
    """
    input:
        CSV = rules.mykrobe_run.output.CSV,
        INFORMS_mykrobe = rules.mykrobe_run.output.INFORMS
    output:
        HTML = 'mykrobe/report/html.iob'# mykrobe.OUTPUT_REPORT
    params:
        show_amr = config['mykrobe'].get('show_amr', True),
        title = config['mykrobe'].get('title', 'Lineage information')
    run:
        from camel.app.tools.mykrobe.mykrobereporter import MykrobeReporter

        reporter = MykrobeReporter()
        if params.show_amr is False:
            reporter.update_parameters(
                show_amr=False,
                custom_header = params.title
            )
        else:
            reporter.update_parameters(custom_header = params.title)

        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        snakemakeutils.add_io_inputs(reporter, input)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule mykrobe_report_empty:
    """
    Creates an empty Mykrobe report when the analysis is disabled.
    """
    output:
        HTML = 'mykrobe/report/html-empty.iob' # mykrobe.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Mykrobe', Path(output.HTML), 3)

rule mykrobe_create_summary:
    """
    Creates the tabular summary output for the Mykrobe assay.
    """
    input:
        INFORMS_mykrobe = rules.mykrobe_run.output.INFORMS
    output:
        TSV = 'mykrobe/summary_mykrobe.{ext}' # mykrobe.OUTPUT_SUMMARY
    params:
        show_amr = config['mykrobe'].get('show_amr', True),
        ext = lambda wildcards: wildcards.ext
    run:
        # Collect informs
        informs = snakemakeutils.load_object(Path(input.INFORMS_mykrobe))

        # Create TSV output
        data_summary = [
            ('mykrobe_phylo_group', informs['phylo_group']),
            ('mykrobe_species', informs['species']),
            ('mykrobe_lineage', informs['lineage']),
        ]
        if params.show_amr is True:
            data_summary.append(('mykrobe_drug_susceptibility', informs['drug_susceptibility']))
        data_summary.extend([
            ('mykrobe_tool_version', informs['_name_full']),
            ('mykrobe_db_version', informs['db_version']),
        ])
        snakemakeutils.export_summary(data_summary, Path(output.TSV), str(params.ext), 'mykrobe')
