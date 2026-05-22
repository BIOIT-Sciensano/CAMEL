import logging
from pathlib import Path

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.viralconsensuspipeline.snakefile import preprocess, nextclade3


rule nextclade3_detect_subtype_mash:
    """
    Uses mash to determine the best matching subtype.
    """
    input:
        FASTQ = preprocess.OUTPUT_FASTQ if config['input']['type'] != 'fasta' else [],
        FASTA = 'fasta.io' if config['input']['type'] == 'fasta' else []
    output:
        TSV = 'nextclade/subtype_determination/mash/tsv.io',
        INFORMS = 'nextclade/subtype_determination/mash/informs.io'
    params:
        dir_ = 'nextclade/subtype_determination/mash',
        input_type = config['input']['type'],
        db = config['nextclade'].get('db_mash')
    run:
        from camel.app.core.errors import ToolExecutionError
        from camelcore.app.io.tooliofile import ToolIOFile
        from camel.app.tools.mash.mashscreen import MashScreen
        from camel.app.scriptutils.basepipe.fastqinput import FastqInput

        mash_screen = MashScreen()
        if params.input_type != 'fasta':
            fq_in = FastqInput.from_fq_dict(Path(input.FASTQ), params.input_type)
            fq_all = []
            for _, tool_io_files in fq_in.to_fq_dict().items():
                fq_all.extend([ToolIOFile(io.path) for io in tool_io_files])
            mash_screen.add_input_files({'FASTQ': fq_all})
        else:
            mash_screen.add_input_files({'FASTA': snakemakeutils.load_object(Path(input.FASTA))})
        path_db = next(Path(params.db).glob('*.msh'))
        logging.info(f'Mash database found: {path_db}')
        mash_screen.add_input_files({'DB': [ToolIOFile(path_db)]})
        try:
            step = Step(rule_name=str(rule), tool=mash_screen, dir_=Path(Path(params.dir_)))
            step.run()
            snakemakeutils.dump_io_outputs(mash_screen, output)
        except ToolExecutionError as err:
            logging.info(f'Error executing {mash_screen.name}')
            snakemakeutils.dump_object([], Path(output.TSV))
            snakemakeutils.dump_object(mash_screen.informs, Path(output.INFORMS))

checkpoint nextclade3_detect_subtype_report:
    """
    Creates the output report with the mash results.
    """
    input:
        TSV = rules.nextclade3_detect_subtype_mash.output.TSV if config['nextclade'].get('dbs') is None else [],
        INFORMS_mash = rules.nextclade3_detect_subtype_mash.output.INFORMS if config['nextclade'].get('dbs') is None else []
    output:
        HTML = 'nextclade/subtype_determination/report/html.iob',
        INFORMS = 'nextclade/subtype_determination/report/informs.io'
    params:
        dir_ = 'nextclade/subtype_determination/report',
        db = config['nextclade'].get('db_mash')
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.nextclade3.nextcladesubtypereporter import NextcladeSubTypeReporter
        from camel.app.core.snakemake import snakepipelineutils

        if params.db is not None:
            reporter = NextcladeSubTypeReporter()
            snakemakeutils.add_io_inputs(reporter, input)
            reporter.add_input_files({'DB': [ToolIODirectory(Path(params.db))]})
            reporter.run(Path(params.dir_))
            snakemakeutils.dump_io_outputs(reporter, output)
        else:
            snakepipelineutils.create_empty_report_section('Subtype determination', Path(output.HTML))
            snakemakeutils.dump_object({}, Path(output.INFORMS))

rule nextclade3_detect_subtype_report_empty:
    """
    Creates an empty output report when subtype detection is not performed.
    """
    output:
        HTML = 'nextclade/subtype_determination/report/html-empty.iob'
    params:
        dir_ = 'nextclade/subtype_determination/report'
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Subtype determination', Path(output.HTML))

rule nextclade3_extract_segment:
    """
    Extracts the given segment from the consensus FASTA file.
    """
    input:
        FASTA = 'fasta.io'
    output:
        FASTA = 'nextclade/{segment}/input/fasta.io'
    params:
        dir_ = lambda wildcards: f'nextclade/{wildcards.segment}/input',
        segment_name = lambda wildcards: wildcards.segment,
        input_type = config['input']['type']
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from camelcore.app.io.tooliofile import ToolIOFile

        # Retrieve the sequence of the corresponding segment
        fasta_in = snakemakeutils.load_object(Path(input.FASTA))[0].path
        with fasta_in.open() as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
            if len(seqs) == 1:
                seq_segment = seqs[0]
            else:
                try:
                    seq_segment = next(
                        s for s in seqs if s.id.lower().endswith(f'-{params.segment_name}'.lower()))
                except StopIteration:
                    segments = [s.id.split('-')[-1].lower() for s in seqs]
                    logger.warning(f"Cannot find segment: {params.segment_name} (found: {', '.join(segments)}), using empty sequence")
                    seq_segment = SeqIO.SeqRecord(seq=Seq(''), id=str(params.segment_name), description="")

        # Save the output file
        path_fasta_out = Path(str(params.dir_), fasta_in.name.replace('.fasta', f'-{params.segment_name}.fasta'))
        with path_fasta_out.open('w') as handle:
            SeqIO.write([seq_segment], handle, 'fasta')
        snakemakeutils.dump_object([ToolIOFile(path_fasta_out)], Path(output.FASTA))

rule nextclade3_run:
    """
    Runs Nextclade on a single sample.
    """
    input:
        FASTA = rules.nextclade3_extract_segment.output.FASTA,
        DB = lambda wildcards: nextclade3.get_nextclade_db(wildcards, checkpoints, segment=wildcards.segment, config=config)
    output:
        TSV = 'nextclade/{segment}/tsv.io',
        INFORMS = 'nextclade/{segment}/informs.io'
    params:
        dir_ = lambda wildcards: f'nextclade/{wildcards.segment}'
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.nextclade3.nextclade3 import Nextclade3

        # Check if database input is valid
        # db_in = snakemakeutils.load_object(Path(input.DB))
        db_in = [ToolIODirectory(Path(str(input.DB)))]
        if len(db_in) == 0:
            snakemakeutils.dump_object([], Path(output.TSV))
            snakemakeutils.dump_object([], Path(output.INFORMS))
        else:
            # Run nextclade
            nextclade = Nextclade3()
            nextclade.add_input_files({
                'FASTA': snakemakeutils.load_object(Path(input.FASTA)),
                'DB': [ToolIODirectory(Path(str(input.DB)))]
            })
            step = Step(rule_name=str(rule), tool=nextclade, dir_=Path(str(params.dir_)))
            step.run()
            snakemakeutils.dump_io_outputs(nextclade, output)

rule nextclade3_reporter:
    """
    Creates an output report for nextclade.
    """
    input:
        TSV = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, key='TSV', config=config),
        INFORMS_nextclade = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, key='INFORMS', config=config)
    output:
        HTML = 'nextclade/report/html.iob' # nextclade3.OUTPUT_REPORT
    params:
        dir_ = 'nextclade/report',
        name = config['input']['sample_name'],
        capitalize_segment_names = config['nextclade'].get('capitalize', False)
    run:
        from camel.app.tools.nextclade3.nextclade3reporter import Nextclade3Reporter

        reporter = Nextclade3Reporter()

        # TSV input
        inputs_tsv = []
        # noinspection PyTypeChecker
        for path_io in [Path(x) for x in input.TSV]:
            inputs_tsv.append(snakemakeutils.load_object(path_io)[0])

        # Nextclade informs
        informs_nextclade = []
        # noinspection PyTypeChecker
        for path_io in [Path(x) for x in input.INFORMS_nextclade]:
            informs_nextclade.append(snakemakeutils.load_object(path_io))

        # Run the reporter
        reporter.add_input_files({'TSV': inputs_tsv})
        reporter.add_input_informs({'nextclade': informs_nextclade})
        reporter.update_parameters(name=params.name, capitalize_segment_names=params.capitalize_segment_names)
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule nextclade_dump_informs:
    """
    Dumps the nextclade informs for the first segment so only a single command is reported, even when multiple segments
    are analyzed.
    """
    input:
        INFORMS_nextclade = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, key='INFORMS', config=config)
    output:
        INFORMS = 'nextclade/informs.io'
    run:
        # noinspection PyUnresolvedReferences
        try:
            informs_first = snakemakeutils.load_object(Path(input.INFORMS_nextclade[0]))
        except IndexError:
            informs_first = []
        snakemakeutils.dump_object(informs_first, Path(output.INFORMS))

rule nextclade3_report_empty:
    """
    Creates an empty report when nextclade is disabled.
    """
    output:
        HTML = 'nextclade/report/html-empty.iob' # nextclade3.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        snakepipelineutils.create_empty_report_section('Nextclade', Path(output.HTML))

rule nextclade3_create_summary:
    """
    Creates the summary output for the nextclade workflow.
    """
    input:
        TSV_nextclade = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, 'TSV', config),
        INFORMS_nextclade = rules.nextclade_dump_informs.output.INFORMS,
        INFORMS_nextclade_all = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, 'INFORMS', config),
        INFORMS_subtype_determination = 'nextclade/subtype_determination/report/informs.io' if (config['nextclade'].get('db') is None) and ('nextclade' in config['analyses_selected']) else []
    output:
        FILE = 'nextclade/summary_nextclade.{ext}' # nextclade3.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        import pandas as pd

        # create the results dictionary:
        assay = 'nextclade'
        data_summary = []

        # Add tool version
        # noinspection PyUnresolvedReferences
        informs_first = snakemakeutils.load_object(Path(input.INFORMS_nextclade))
        data_summary.append((f'{assay}_version', informs_first['_version']))

        # Add subtype determination informs (if available)
        if input.INFORMS_subtype_determination:
           informs_subtype_determination = snakemakeutils.load_object(Path(input.INFORMS_subtype_determination))
           data_summary.append((f'{assay}_detected_subtype', informs_subtype_determination.get('subtype')))

        # Nextclade output
        keys_kept = [
            'qc.frameShifts.status',
            'qc.mixedSites.status',
            'qc.mixedSites.status',
            'qc.overallScore',
            'qc.overallStatus',
            'qc.privateMutations.status',
            'qc.stopCodons.status',
        ]
        # noinspection PyTypeChecker
        for path_io, path_informs in zip(
                [Path(x) for x in input.TSV_nextclade],
                [Path(x) for x in input.INFORMS_nextclade_all]):

            # Determine prefix
            path_tsv = snakemakeutils.load_object(path_io)[0].path
            segment_name = path_tsv.parents[1].name
            # noinspection PyTypeChecker
            base_key = assay if len(input.TSV_nextclade) == 1 else f'{assay}_{segment_name}'

            # Add results data
            data_in = pd.read_table(path_tsv).fillna('-')
            for key in keys_kept:
                data_summary.append((
                    f'{base_key}_{key}',
                    data_in.iloc[0][key]
                ))

            # Add informs data
            informs = snakemakeutils.load_object(path_informs)
            data_summary.append((f'{base_key}_reference', informs['db']['reference']))
            data_summary.append((f'{base_key}_version', informs['db']['version']))

            # Add metadata
            if informs['db'].get('metadata_columns') is not None:
                for row in informs['db']['metadata_columns']:
                    data_summary.append((f"{base_key}_{row['key']}", data_in.iloc[0][row['key']]))

        # Save in TSV format
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'nextclade')
