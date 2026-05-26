import shutil
from pathlib import Path
from typing import Any

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.snakefiles import sequence_typing


def aggregate_input(wc: Any, cp: Any, conf: dict, get_informs: bool) -> list[str]:
    """
    Determines which output files need to be generated depending on the scheme and detection method.
    :param wc: Snakemake wildcards
    :param cp: Snakemake checkpoints
    :param conf: Snakemake config
    :param get_informs: If True, the informs are returned
    :return: List of output files
    """
    # Parse scheme informs
    path_metadata = cp.typing_extract_schema_info.get(**wc).output.INFORMS
    informs_ = snakemakeutils.load_object(Path(path_metadata))

    # Fallback when database is missing
    if 'loci' not in informs_:
        return []

    # Determine which loci are present
    has_dna_loci = any(info['type'] == 'DNA' for _, info in informs_['loci'].metadata_by_locus_name.items())
    has_peptide_loci = any(info['type'] == 'peptide' for _, info in informs_['loci'].metadata_by_locus_name.items())

    # Create output files
    basename = 'hits.iob' if not get_informs else 'informs.io'
    output_files = []
    if has_dna_loci:
        method = sequence_typing.get_detection_method(conf, wc.scheme, 'DNA')
        if method == 'mist':
            output_files.append(f'typing/{wc.scheme}/mist/DNA/{basename}')
        else:
            output_files.append(f'typing/{wc.scheme}/type_async/DNA/{method}/{basename}')
    if has_peptide_loci:
        output_files.append(f'typing/{wc.scheme}/type_async/peptide/blast/{basename}')
    return [str(x) for x in output_files]

checkpoint typing_extract_schema_info:
    """
    Extracts the metadata for a scheme.
    """
    output:
        INFORMS = 'typing/{scheme}/scheme_info/informs.iob', # sequence_typing.OUTPUT_DB_INFORMS
        TSV = 'typing/{scheme}/scheme_info/tsv-profiles.io'
    params:
        db = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme]['path']
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.sequence_typing.typingdbloader import TypingDBLoader
        db_loader = TypingDBLoader()
        path_db = Path(str(params.db))
        if path_db.exists():
            db_loader.add_input_files({'DIR': [ToolIODirectory(Path(str(params.db)))]})
        else:
            db_loader.update_parameters(db_name=path_db.name)
        step = Step(rule_name=str(rule), tool=db_loader, dir_=snakemakeutils.get_rule_dir(output), wildcards=wildcards)
        step.run()
        snakemakeutils.dump_object(db_loader.informs, Path(output.INFORMS))
        snakemakeutils.dump_object(
            db_loader.tool_outputs['TSV'] if 'TSV' in db_loader.tool_outputs else [], Path(output.TSV))

rule typing_async:
    """
    Types all loci asynchronously.
    Supported detection methods: blast, KMA.
    """
    input:
        FASTA = lambda wildcards: str(sequence_typing.INPUT_FASTA) if wildcards.detection_method == 'blast' else [],
        FASTQ = lambda wildcards: 'fq_dict.io' if wildcards.detection_method in ('kma',) else [],
        DIR = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme]['path']
    output:
        IO = 'typing/{scheme}/type_async/{locus_type}/{detection_method}/hits.iob',
        INFORMS = 'typing/{scheme}/type_async/{locus_type}/{detection_method}/informs.io'
    params:
        detection_method = lambda wildcards: wildcards.detection_method,
        locus_type = lambda wildcards: wildcards.locus_type,
        blastn_task = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme].get('blastn_task', 'megablast'),
        read_type = 'PE' if config['input'].get('type', 'illumina') == 'illumina' else 'SE'
    threads: 8
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.pipelines.sequence_typing.typeasync import TypeAsync
        from camelcore.app.io.tooliodirectory import ToolIODirectory

        # Create the working directory
        dir_working = snakemakeutils.get_rule_dir(output)
        dir_working.mkdir(parents=True, exist_ok=True)

        # Run the tool
        typer = TypeAsync()
        step = Step(rule_name=str(rule), tool=typer, dir_=dir_working, wildcards=wildcards)
        if params.detection_method == 'blast':
            typer.add_input_files({'FASTA': snakemakeutils.load_object(Path(str(input.FASTA)))})
            typer.update_parameters(blastn_task=str(params.blastn_task))
        elif params.detection_method == 'kma':
            fastq_input = snakepipelineutils.extract_fq_input(
                Path(str(input.FASTQ)), key_pe='FASTQ_PE', key_se='FASTQ_SE', read_type=params.read_type)
            typer.add_input_files(fastq_input)
        else:
            raise ValueError(f'Invalid detection method: {params.detection_method}')
        typer.add_input_files({'DIR': [ToolIODirectory(Path(str(input.DIR)))]})
        typer.update_parameters(
            detection_method=str(params.detection_method), locus_type=str(params.locus_type), threads=threads)
        step.run()

        # Save the tool output
        snakemakeutils.dump_object(typer.tool_outputs['VAL_hits'], Path(output.IO))
        snakemakeutils.dump_object(typer.informs, Path(output.INFORMS))

rule typing_mist:
    """
    Performs sequence typing using MiST.
    """
    input:
        FASTA = sequence_typing.INPUT_FASTA,
        DIR = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme]['path']
    output:
        JSON = 'typing/{scheme}/mist/DNA/json.io',
        INFORMS = 'typing/{scheme}/mist/DNA/informs.io'
    threads: 4
    run:
        from camelcore.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.mist.mistcall import MiSTCall
        mist_call = MiSTCall()
        mist_call.add_input_files({
            'FASTA': snakemakeutils.load_object(Path(input.FASTA)),
            'DB': [ToolIODirectory(Path(str(input.DIR), 'mist'))]
        })
        mist_call.update_parameters(threads=threads)
        step = Step(str(rule), mist_call, dir_=snakemakeutils.get_rule_dir(output))
        step.run()
        snakemakeutils.dump_io_outputs(mist_call, output)

rule typing_mist_extract_hits:
    """
    Extracts sequence typing hits from the MiST output.
    """
    input:
        JSON = rules.typing_mist.output.JSON
    output:
        IO = 'typing/{scheme}/mist/DNA/hits.iob'
    run:
        from camelcore.app.io.tooliovalue import ToolIOValue
        from camel.app.toolkits.sequencetyping.typingmisthit import TypingMiSTHit
        path_json = snakemakeutils.load_object(Path(input.JSON))[0].path
        hits = TypingMiSTHit.parse_mist_json(path_json)
        snakemakeutils.dump_object([ToolIOValue(h) for h in hits], Path(output.IO))

rule typing_collect_hits:
    """
    Helper rule to collect the hits from the typing_async rule.
    """
    input:
        lambda wildcards: aggregate_input(wildcards, checkpoints, config, get_informs=False)
    output:
        IO = 'typing/{scheme}/hits/{locus_type}/hits.iob'
    params:
        locus_type = lambda wildcards: wildcards.locus_type
    run:
        # Retrieve the hits of the matching locus type
        for path_in in [Path(x) for x in input]:
            for segment in path_in.parents:
                if segment.name == params.locus_type:
                    shutil.copyfile(path_in, output.IO)
                    return

        # No matching hits found -> empty output
        snakemakeutils.dump_object([], Path(output.IO))

rule typing_add_allele_page_url:
    """
    Adds the locus URL to the typing hits (if available).
    """
    input:
        IO = rules.typing_collect_hits.output.IO,
        INFORMS_scheme = rules.typing_extract_schema_info.output.INFORMS
    output:
        IO = 'typing/{scheme}/hits/{locus_type}/hits-url.iob'
    run:
        from camelcore.app.utils import fileutils

        # Load the informs
        metadata_by_locus_name = snakemakeutils.load_object(
            Path(str(input.INFORMS_scheme)))['loci'].metadata_by_locus_name

        # Add the allele url to the hit
        hits = snakemakeutils.load_object(Path(input.IO))
        for hit in hits:
            locus_key = fileutils.make_valid(hit.value.locus).lower()
            locus_metadata = metadata_by_locus_name[locus_key]
            hit.value.set_allele_page_url_template(locus_metadata.get('allele_page_url'))

        # Export hits
        snakemakeutils.dump_object(hits, Path(output.IO))

rule typing_export_tsv:
    """
    Exports the detected typing hits in TSV format.
    """
    input:
        HITS = rules.typing_collect_hits.output.IO
    output:
        TSV = 'typing/{scheme}/export_tsv/{locus_type}/tsv.io' # sequence_typing.OUTPUT_TSV
    params:
        output_filename = lambda wildcards: f"typing-{wildcards.scheme}-{wildcards.locus_type}-{config['input']['sample_name']}.tsv"
    run:
        from camel.app.tools.pipelines.sequence_typing.exporttsv import ExportTSV
        export_tsv = ExportTSV()
        snakemakeutils.add_io_inputs(export_tsv, input)
        step = Step(str(rule), export_tsv, dir_=snakemakeutils.get_rule_dir(output))
        export_tsv.update_parameters(output_filename=str(params.output_filename))
        step.run()
        snakemakeutils.dump_object(
            export_tsv.tool_outputs['TSV'] if 'TSV' in export_tsv.tool_outputs else [], Path(output.TSV))

rule typing_dump_commands:
    """
    Collects the commands used to perform the typing.
    """
    input:
        INFORMS_detection = lambda wildcards: aggregate_input(wildcards, checkpoints, config, get_informs=True),
        INFORMS_scheme = rules.typing_extract_schema_info.output.INFORMS
    output:
        INFORMS = 'typing/{scheme}/commands/informs.io' # sequence_typing.OUTPUT_INFORMS
    params:
        scheme = lambda wildcards: wildcards.scheme
    run:
        informs_scheme = snakemakeutils.load_object(Path(input.INFORMS_scheme))
        informs_commands = []
        # noinspection PyTypeChecker
        for path_informs in [Path(x) for x in input.INFORMS_detection]:
            data = snakemakeutils.load_object(path_informs)
            if len(data) == 0 or '_command' not in data:
                continue
            data['_tag'] = f"Typing - {informs_scheme['title']}"
            informs_commands.append(data)
        snakemakeutils.dump_object(informs_commands, Path(output.INFORMS))

rule typing_detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
        hits_nucl = lambda wildcards: str(rules.typing_collect_hits.output.IO).format(scheme=wildcards.scheme, locus_type='DNA'),
        hits_pept = lambda wildcards: str(rules.typing_collect_hits.output.IO).format(scheme=wildcards.scheme, locus_type='peptide'),
        TSV = rules.typing_extract_schema_info.output.TSV
    output:
        TSV = 'typing/{scheme}/detect_st/tsv.io',
        INFORMS = 'typing/{scheme}/detect_st/informs-st.io'
    params:
        write_all_matches = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme].get('write_all_matches', False)
    run:
        from camel.app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        data_profiles = snakemakeutils.load_object(Path(input.TSV))
        if len(data_profiles) == 0:
            snakemakeutils.dump_object([], Path(output.INFORMS))
            snakemakeutils.dump_object([], Path(output.TSV))
        else:
            sequence_type_detector = SequenceTypeDetector()
            snakemakeutils.add_io_inputs(sequence_type_detector, input)
            step = Step(rule_name=str(rule), tool=sequence_type_detector, dir_=snakemakeutils.get_rule_dir(output), wildcards=wildcards)
            sequence_type_detector.update_parameters(allele_wildcards='N', allele_absent_symbol='0')
            if params.write_all_matches:
                 sequence_type_detector.update_parameters(write_tsv='True', output_filename='profile_matches.tsv')
            step.run()
            snakemakeutils.dump_object(sequence_type_detector.informs, Path(output.INFORMS))
            snakemakeutils.dump_object(
                sequence_type_detector.tool_outputs['TSV'] if params.write_all_matches else [], Path(output.TSV))

rule typing_get_cgmlst_stats:
    """
    Retrieves the cgMLST stats for the given scheme.
    Only DNA loci are considered.
    """
    input:
        HITS = lambda wildcards: str(rules.typing_collect_hits.output.IO).format(
            scheme=wildcards.scheme, locus_type='DNA')
    output:
        INFORMS = 'typing/{scheme}/stats/informs.io'
    params:
        scheme_name = lambda wildcards: wildcards.scheme
    run:
        all_hits = snakemakeutils.load_object(Path(str(input.HITS)))
        nb_found = len([v for v in all_hits if v.value.is_perfect_hit() or v.value.is_new_allele()])
        snakemakeutils.dump_object({
            'hits_found': nb_found,
            'nb_of_loci': len(all_hits),
            'scheme_name': params.scheme_name},
            Path(output.INFORMS))

# noinspection PyUnresolvedReferences
rule typing_create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV_nucl = rules.typing_export_tsv.output.TSV.format(locus_type='DNA', scheme='{scheme}'),
        TSV_pept = rules.typing_export_tsv.output.TSV.format(locus_type='peptide', scheme='{scheme}'),
        INFORMS_scheme = rules.typing_extract_schema_info.output.INFORMS,
        INFORMS_ST = rules.typing_detect_sequence_type.output.INFORMS,
        hits_nucl = rules.typing_add_allele_page_url.output.IO.format(locus_type='DNA', scheme='{scheme}'),
        hits_pept = rules.typing_add_allele_page_url.output.IO.format(locus_type='peptide', scheme='{scheme}')
    output:
        VAL_HTML = 'typing/{scheme}/report/html.iob' # sequence_typing.OUTPUT_REPORT
    params:
        sample_name = config['input']['sample_name'],
        detection_method = lambda wildcards: sequence_typing.get_detection_method(config, wildcards.scheme, 'DNA'),
        config_data = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme]
    run:
        from camelcore.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping()

        # Add ST informs (if not empty)
        data_informs = snakemakeutils.load_object(Path(input.INFORMS_ST))
        if len(data_informs) > 0:
            reporter.add_input_informs({'ST': snakemakeutils.load_object(Path(str(input.INFORMS_ST)))})

        # Add other inputs and parameters
        snakemakeutils.add_io_inputs(reporter, input, excluded_keys=['INFORMS_ST'])
        reporter.add_input_files({'VAL_SAMPLE': [ToolIOValue(params.sample_name)]})
        if params.detection_method != config['sequence_typing']['options']['method']:
            reporter.update_parameters(forced_detection_method=str(params.detection_method))
        # noinspection PyUnresolvedReferences
        if params.config_data.get('hidden', False) is True:
            reporter.update_parameters(hidden=True)

        # Optional message
        if 'message' in params.config_data:
            reporter.update_parameters(message=params.config_data['message']['content'])
        if ('message' in params.config_data) and ('category' in params.config_data['message']):
            reporter.update_parameters(message_category=params.config_data['message']['category'])

        # Run the reporter
        step = Step(rule_name=str(rule), tool=reporter, dir_=snakemakeutils.get_rule_dir(output), wildcards=wildcards)
        step.run()
        snakemakeutils.dump_io_outputs(reporter, output)

rule typing_create_report_empty:
    """
    Creates an empty sequence typing report when the analysis is disabled.
    """
    input:
        INFORMS = rules.typing_extract_schema_info.output.INFORMS
    output:
        HTML = 'typing/{scheme}/report/html-empty.iob' # sequence_typing.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.core.snakemake import snakepipelineutils
        informs = snakemakeutils.load_object(Path(str(input.INFORMS)))
        snakepipelineutils.create_empty_report_section(informs['title'], Path(output.HTML))

rule typing_dump_summary_info:
    """
    Dumps the summary information in tabular format.
    """
    input:
        HITS_NUCL = lambda wildcards: str(rules.typing_add_allele_page_url.output.IO).format(scheme=wildcards.scheme, locus_type='DNA'),
        HITS_PEPT = lambda wildcards: str(rules.typing_add_allele_page_url.output.IO).format(scheme=wildcards.scheme, locus_type='peptide'),
        INFORMS_ST = rules.typing_detect_sequence_type.output.INFORMS
    output:
        FILE = 'typing/{scheme}/summary/summary_out.{ext}' # sequence_typing.OUTPUT_SUMMARY
    params:
        scheme_name = lambda wildcards: wildcards.scheme,
        ext = lambda wildcards: wildcards.ext,
        include_hashing = lambda wildcards: config['sequence_typing']['dbs'][wildcards.scheme].get('include_hashing')
    run:
        summary_data = []

        ##########################
        # Sequence type metadata #
        ##########################
        data_st = snakemakeutils.load_object(Path(input.INFORMS_ST))
        if len(data_st) == 0:
            st_metadata = []
        else:
            st_metadata = data_st['metadata']
            st_metadata.append(('percent_detected', data_st['percent_detected']))
        for k, v in st_metadata:
            summary_data.append((f'{params.scheme_name}-{k}', v))

        ########
        # Loci #
        ########
        hits = snakemakeutils.load_object(Path(str(input.HITS_NUCL))) + snakemakeutils.load_object(Path(str(input.HITS_PEPT)))

        if params.ext == 'tsv':
            # TSV format: Each line represents a locus, with table values separated by commas
            for hit in [io.value for io in hits]:
                key = f'{params.scheme_name}-{hit.locus}'
                data_hit = ','.join(hit.to_table_row())
                summary_data.append((key, data_hit))
        elif params.ext == 'json':
            # JSON format: A list of loci with corresponding information as a dictionary
            summary_data.append(('loci', [hit.to_dict(params.include_hashing) for hit in [io.value for io in hits]]))
        else:
            raise ValueError(f'Invalid extension: {params.ext}')

        # Create output
        snakemakeutils.export_summary(summary_data, Path(output.FILE), str(params.ext), str(params.scheme_name))
