import logging
from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.app.tools.mash.mashscreen import MashScreen
from camel.scripts.viralconsensuspipeline.snakefile import refselection


rule ref_selection_mash_screen:
    """
    Calculates similarity to genomes in the reference database using mash screen.
    """
    input:
        IO_FASTQ = 'fq_dict.io',
        DB = lambda wildcards: Path(config['ref_selection']['db'], 'mash', f'{wildcards.segment}.msh') if config['ref_selection']['db'] is not None else []
    output:
        TSV = 'ref_selection/mash_screen/{segment}/tsv.io',
        INFORMS = 'ref_selection/mash_screen/{segment}/informs.io'
    threads: 8
    params:
        input_type = config['input']['type'],
        segment = lambda wildcards: wildcards.segment
    run:
        from camelcore.app.io.tooliofile import ToolIOFile

        # Load FASTQ data
        if params.input_type == 'ont':
            fq_in = [f.path for f in snakemakeutils.load_object(Path(input.IO_FASTQ))['SE']]
        else:
            fq_in = [f.path for f in snakemakeutils.load_object(Path(input.IO_FASTQ))['PE']]

        mash = MashScreen()
        mash.add_input_files({
            'FASTQ': [ToolIOFile(p) for p in fq_in],
            'DB': [ToolIOFile(Path(str(input.DB)))]
        })
        mash.update_parameters(threads=threads)
        step = Step(rule_name=str(rule), tool=mash, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(mash, output)

rule ref_selection_create_fasta:
    """
    Creates a FASTA file by parsing the mash output and selection the best matching reference for each segment.
    """
    input:
        TSV = ([str(rules.ref_selection_mash_screen.output.TSV).format(segment=seg) for seg in refselection.get_segments(
            Path(config['ref_selection']['db']))]) if config['ref_selection']['db'] is not None else [],
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        TSV = 'ref_selection/create_fasta/tsv.io',
        JSON = 'ref_selection/create_fasta/json.io',
        FASTA = 'ref_selection/create_fasta/fasta.io' # refselection.OUTPUT_FASTA
    run:
        import itertools
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.viral_consensus.refselectioncreatefasta import RefSelection
        ref_selection = RefSelection()
        tsv_in = list(itertools.chain(*[snakemakeutils.load_object(Path(io)) for io in input.TSV]))
        ref_selection.add_input_files({'TSV': tsv_in, 'DB': [ToolIODirectory(Path(input.DB))]})
        step = Step(rule_name=str(rule), tool=ref_selection, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(ref_selection, output)

rule ref_selection_report:
    """
    Creates the HTML output report for the reference selection.
    """
    input:
        TSV = rules.ref_selection_create_fasta.output.TSV,
        JSON = rules.ref_selection_create_fasta.output.JSON,
        FASTA = rules.ref_selection_create_fasta.output.FASTA,
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        VAL_HTML = 'ref_selection/report/html.iob' # refselection.OUTPUT_REPORT
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.viral_consensus.reporterrefselection import ReporterRefSelection
        reporter = ReporterRefSelection()
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output))
        snakemakeutils.add_io_inputs(reporter, input, excluded_keys=['DB'])
        reporter.add_input_files({'DB': [ToolIODirectory(Path(input.DB))]})
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule ref_selection_report_empty:
    """
    Creates an empty report when reference genome selection is disabled. 
    """
    output:
        VAL_HTML = 'ref_selection/report/html-empty.iob' # refselection.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Reference selection', Path(output.VAL_HTML))

rule ref_selection_dump_summary_info:
    """
    Creates the summary information for the reference selection.
    """
    input:
        JSON = rules.ref_selection_create_fasta.output.JSON,
        DB = Path(config['ref_selection']['db']) if config['ref_selection'].get('db') is not None else []
    output:
        FILE = 'ref_selection/summary/summary.{ext}' # refselection.OUTPUT_REF_SELECTION_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import json

        # Load JSON input
        path_json = snakemakeutils.load_object(Path(input.JSON))[0].path
        logging.info(f'Reading: {path_json}')
        with path_json.open() as handle:
            data_ref_selection = json.load(handle)

        # Output the selected references
        data_summary = []
        data_summary.append(('ref_selection_database', Path(input.DB).name if input.DB else '-'))
        for segment, ref in data_ref_selection.items():
            data_summary.append((f"ref_selection-{segment}", ref if ref is not None else '-'))
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'ref_selection')
