from pathlib import Path

from snakemake.checkpoints import Checkpoints

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import sequence_typing



def aggregate_input(wc: dict, cp: Checkpoints, conf: dict) -> list[str]:
    """
    Aggregates the input to obtain the required output files.
    :param wc: Snakemake wildcards
    :param cp: Snakemake checkpoints
    :param conf: Snakemake config
    """
    # Parse scheme informs
    path_metadata = cp.typing_extract_schema_info.get(**wc).output.INFORMS
    informs = snakemakeutils.load_object(Path(path_metadata))

    # Determine which loci are present
    has_dna_loci = any(info['type'] == 'DNA' for _, info in  informs['loci'].metadata_by_locus_name.items())
    has_peptide_loci = any(info['type'] == 'peptide' for _, info in  informs['loci'].metadata_by_locus_name.items())

    # Create output files
    output_files = []
    if has_dna_loci:
        method = SequenceTypingUtils.get_detection_method(conf, wc.scheme)
        output_files.append(f'typing/{wc.scheme}/DNA/{method}/hits.iob')
        output_files.append(f'typing/{wc.scheme}/DNA/{method}/informs.io')
    if has_peptide_loci:
        output_files.append(f'typing/{wc.scheme}/peptide/blast/hits.iob')
        output_files.append(f'typing/{wc.scheme}/peptide/blast/informs.io')
    return [str(x) for x in output_files]

checkpoint typing_extract_schema_info:
    """
    Extracts the metadata for a scheme.
    """
    input:
        DIR = lambda wildcards: config['sequence_typing'][wildcards.scheme]['path']
    output:
        INFORMS = 'typing/{scheme}/informs-locus_set.iob'
    params:
        running_dir = lambda wildcards: f'typing/{wildcards.scheme}'
    run:
        from camel.app.io.tooliodirectory import ToolIODirectory
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        locus_set_manager = LocusSetManager()
        locus_set_manager.add_input_files({'DIR': [ToolIODirectory(Path(str(input.DIR)))]})
        step = Step(rule_name=str(rule), tool=locus_set_manager, dir_=Path(str(params.running_dir)), wildcards=wildcards)
        step.run()
        snakemakeutils.dump_object(locus_set_manager.informs, Path(output.INFORMS))

rule typing_pickle_profiles:
    """
    Retrieves the sequence type definitions and converts them to CAMEL IO pickles.
    """
    input:
        DIR = lambda wildcards: config['sequence_typing'][wildcards.scheme]['path']
    output:
        TSV = 'typing/{scheme}/tsv-profiles.io'
    run:
        from camel.app.io.tooliofile import ToolIOFile
        path_profiles = Path(str(input.DIR)) / 'profiles.tsv'
        snakemakeutils.dump_object([ToolIOFile(Path(
            path_profiles))] if path_profiles.exists() else [], Path(output.TSV))

rule typing_async:
    """
    Types all loci asynchronously.
    """
    input:
        FASTA = lambda wildcards: str(sequence_typing.INPUT_FASTA) if wildcards.detection_method == 'blast' else [],
        FASTQ = lambda wildcards: (str('fq_dict.io')) if wildcards.detection_method in ('kma',) else [],
        DIR = lambda wildcards: config['sequence_typing'][wildcards.scheme]['path']
    output:
        IO = 'typing/{scheme}/{locus_type}/{detection_method}/hits.iob',
        INFORMS = 'typing/{scheme}/{locus_type}/{detection_method}/informs.io'
    params:
        dir_working = lambda wildcards: f'typing/{wildcards.scheme}/{wildcards.locus_type}',
        detection_method = lambda wildcards: wildcards.detection_method,
        locus_type = lambda wildcards: wildcards.locus_type,
        blastn_task = lambda wildcards: config['sequence_typing'][wildcards.scheme].get('blastn_task', 'megablast'),
        read_type = 'PE' if config.get('input_type', 'illumina') == 'illumina' else 'SE'
    threads: max(16, workflow.cores * 0.75)
    run:
        from camel.app.tools.pipelines.sequence_typing.typeasync import TypeAsync
        from camel.app.io.tooliodirectory import ToolIODirectory

        # Create working directory
        dir_working = Path(str(params.dir_working))
        dir_working.mkdir(parents=True, exist_ok=True)

        # Run the tool
        typer = TypeAsync()
        step = Step(rule_name=str(rule), tool=typer, dir_=dir_working, wildcards=wildcards)
        if params.detection_method == 'blast':
            typer.add_input_files({'FASTA': snakemakeutils.load_object(Path(str(input.FASTA)))})
            typer.update_parameters(blastn_task=str(params.blastn_task))
        elif params.detection_method in ('kma',):
            fastq_input = SnakePipelineUtils.extracts_fq_input(
                Path(str(input.FASTQ)), key_pe='FASTQ_PE', key_se='FASTQ_SE', read_type=params.read_type)
            typer.add_input_files(fastq_input)
        typer.add_input_files({'DIR': [ToolIODirectory(Path(str(input.DIR)))]})
        typer.update_parameters(
            detection_method=str(params.detection_method), locus_type=str(params.locus_type), threads=threads)
        step.run()

        # Save the tool output
        snakemakeutils.dump_object(typer.tool_outputs['VAL_hits'], Path(output.IO))
        snakemakeutils.dump_object(typer.informs, Path(output.INFORMS))

rule typing_get_hits:
    """
    Selects the hits output based on the detection method in the config.
    """
    input:
        lambda wildcards: aggregate_input(wildcards, checkpoints, config)
    output:
        HITS_NUCL = 'typing/{scheme}/DNA/hits.iob',
        HITS_PEPT = 'typing/{scheme}/peptide/hits.iob',
        INFORMS = 'typing/{scheme}/informs_detection.io'
    run:
        informs = []
        hits_by_type = {'DNA': [], 'peptide': []}
        for path_io in [Path(x) for x in input]:
            locus_type = path_io.parents[1].name
            if path_io.stem == 'hits':
                hits_by_type[locus_type].extend(snakemakeutils.load_object(path_io))
            elif path_io.stem == 'informs':
                informs.append(snakemakeutils.load_object(path_io))
            else:
                raise ValueError(f"Invalid input file: {path_io}")
        snakemakeutils.dump_object(hits_by_type['DNA'], Path(output.HITS_NUCL))
        snakemakeutils.dump_object(hits_by_type['peptide'], Path(output.HITS_PEPT))
        snakemakeutils.dump_object(informs, Path(output.INFORMS))

rule typing_dump_commands:
    input:
        INFORMS_detection = rules.typing_get_hits.output.INFORMS,
        INFORMS_scheme = str('typing/{scheme}/informs-locus_set.iob')
    output:
        INFORMS = sequence_typing.OUTPUT_INFORMS
    params:
        scheme = lambda wildcards: wildcards.scheme
    run:
        informs_scheme = snakemakeutils.load_object(Path(input.INFORMS_scheme))
        informs_commands = []
        for informs in snakemakeutils.load_object(Path(input.INFORMS_detection)):
            informs['_tag'] = f"Typing - {informs_scheme['title']}"
            informs_commands.append(informs)
        snakemakeutils.dump_object(informs_commands, Path(output.INFORMS))

rule typing_export_hits_tabular:
    """
    Creates a tabular output for the detected hits.
    """
    input:
        hits = str('typing/{scheme}/{locus_type}/hits.iob')
    output:
        TSV = sequence_typing.OUTPUT_TSV
    params:
        dir_ = lambda wildcards: f'typing/{wildcards.scheme}/{wildcards.locus_type}/tabular',
        sample_name = FileSystemHelper.make_valid(config['sample_name']),
        scheme = lambda wildcards: FileSystemHelper.make_valid(wildcards.scheme),
        locus_type = lambda wildcards: wildcards.locus_type
    run:
        from camel.app.io.tooliofile import ToolIOFile
        import pandas as pd

        hits = [h.value for h in snakemakeutils.load_object(Path(input.hits))]

        # No hits detected -> no TSV file is generated
        if len(hits) == 0:
            snakemakeutils.dump_object([], Path(output.TSV))
            return

        # Export TSV file
        path_out = Path(str(params.dir_), f'typing-{params.scheme}-{params.locus_type}-{params.sample_name}.tsv')
        data_hits = pd.DataFrame(data=[h.to_table_row() for h in hits], columns=hits[0].table_column_names())
        data_hits.to_csv(path_out, sep='\t', index=False)
        snakemakeutils.dump_object([ToolIOFile(path_out)], Path(output.TSV))

        # Export hashed TSV file if there are any novel alleles
        if not any([h.is_new_allele() for h in hits]):
            return
        path_out_hash = Path(str(params.dir_), f'typing-{params.scheme}-{params.locus_type}-{params.sample_name}-hashes.tsv')
        data_hits = pd.DataFrame(
            data=[h.to_table_row(hash_allele_ids=True) for h in hits], columns=hits[0].table_column_names())
        data_hits.to_csv(path_out_hash, sep='\t', index=False)

rule typing_detect_sequence_type:
    """
    Detects the sequence type based on the detected alleles.
    """
    input:
        hits_nucl = rules.typing_get_hits.output.HITS_NUCL,
        hits_pept = rules.typing_get_hits.output.HITS_PEPT,
        TSV = rules.typing_pickle_profiles.output.TSV
    output:
        TSV = 'typing/{scheme}/detect_st/tsv.io',
        INFORMS = 'typing/{scheme}/detect_st/informs-st.io'
    params:
        running_dir = lambda wildcards: f'typing/{wildcards.scheme}/detect_st',
        write_all_matches = lambda wildcards: config['sequence_typing'][wildcards.scheme].get('write_all_matches', False)
    run:
        from camel.app.tools.pipelines.sequence_typing.sequencetypedetector import SequenceTypeDetector
        data_profiles = snakemakeutils.load_object(Path(input.TSV))
        if len(data_profiles) == 0:
            snakemakeutils.dump_object([], Path(output.INFORMS))
            snakemakeutils.dump_object([], Path(output.TSV))
        else:
            sequence_type_detector = SequenceTypeDetector()
            snakemakeutils.add_pickle_inputs(sequence_type_detector, input)
            step = Step(rule_name=str(rule), tool=sequence_type_detector, dir_=Path(str(params.running_dir)), wildcards=wildcards)
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
        HITS = rules.typing_get_hits.output.HITS_NUCL
    output:
        INFORMS = 'typing/{scheme}/stats/informs.io'
    params:
        scheme_name = lambda wildcards: wildcards.scheme
    run:
        all_hits = snakemakeutils.load_object(Path(input.HITS))
        nb_found = len([v for v in all_hits if v.value.is_perfect_hit() or v.value.is_new_allele()])
        snakemakeutils.dump_object(
            {'hits_found': nb_found,
             'nb_of_loci': len(all_hits),
             'scheme_name': params.scheme_name},
            Path(output.INFORMS))

rule typing_add_allele_page_url:
    """
    This steps add the locus url to the detected hits.
    """
    input:
        HITS_NUCL = rules.typing_get_hits.output.HITS_NUCL,
        HITS_PEPT = rules.typing_get_hits.output.HITS_PEPT,
        INFORMS_scheme = lambda wildcards: checkpoints.typing_extract_schema_info.get(scheme=wildcards.scheme).output.INFORMS
    output:
        HITS_NUCL = 'typing/{scheme}/DNA/hits-url.iob',
        HITS_PEPT = 'typing/{scheme}/peptide/hits-url.iob'
    run:
        for key in ('NUCL', 'PEPT'):
            hits = snakemakeutils.load_object(Path(input.get(f'HITS_{key}')))

            # Load the informs
            metadata_by_locus_name = snakemakeutils.load_object(
                Path(str(input.INFORMS_scheme)))['loci'].metadata_by_locus_name

            # Add the allele url to the hit
            for hit in hits:
                locus_key = FileSystemHelper.make_valid(hit.value.locus).lower()
                locus_metadata = metadata_by_locus_name[locus_key]
                hit.value.set_allele_page_url_template(locus_metadata.get('allele_page_url'))

            # Export hits
            snakemakeutils.dump_object(hits, Path(output.get(f'HITS_{key}')))

rule typing_create_report:
    """
    Creates a report with the sequence typing output.
    """
    input:
        TSV_nucl = rules.typing_export_hits_tabular.output.TSV.format(locus_type='DNA', scheme='{scheme}'),
        TSV_pept = rules.typing_export_hits_tabular.output.TSV.format(locus_type='peptide', scheme='{scheme}'),
        INFORMS_scheme = lambda wildcards: checkpoints.typing_extract_schema_info.get(scheme=wildcards.scheme).output.INFORMS,
        INFORMS_ST = rules.typing_detect_sequence_type.output.INFORMS,
        hits_nucl = rules.typing_add_allele_page_url.output.HITS_NUCL,
        hits_pept = rules.typing_add_allele_page_url.output.HITS_PEPT
    output:
        VAL_HTML = 'typing/{scheme}/html.iob' #
    # sequence_typing.OUTPUT_REPORT
    params:
        running_dir = lambda wildcards: f'typing/{wildcards.scheme}',
        sample_name = config['sample_name'],
        detection_method = lambda wildcards: SequenceTypingUtils.get_detection_method(config, wildcards.scheme),
        config_data = lambda wildcards: config['sequence_typing'][wildcards.scheme]
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping
        reporter = HtmlReporterTyping()

        # Add ST informs (if not empty)
        data_informs = snakemakeutils.load_object(Path(input.INFORMS_ST))
        if len(data_informs) > 0:
            reporter.add_input_informs({'ST': snakemakeutils.load_object(Path(str(input.INFORMS_ST)))})

        # Add other inputs and parameters
        snakemakeutils.add_pickle_inputs(reporter, input, excluded_keys=['INFORMS_ST'])
        reporter.add_input_files({'VAL_SAMPLE': [ToolIOValue(params.sample_name)]})
        if params.detection_method != config['detection_method']:
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
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(str(params.running_dir)), wildcards=wildcards)
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule typing_create_report_empty:
    """
    Creates an empty sequence typing report when the analysis is disabled.
    """
    input:
        INFORMS_Scheme = lambda wildcards: checkpoints.typing_extract_schema_info.get(scheme=wildcards.scheme).output.INFORMS
    output:
        VAL_HTML = 'typing/{scheme}/html-empty.iob' # sequence_typing.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.io.tooliovalue import ToolIOValue
        informs = snakemakeutils.load_object(Path(str(input.INFORMS_Scheme)))
        section = HtmlReportSection(informs['title'], 3)
        section.add_paragraph('Analysis disabled')
        snakemakeutils.dump_object([ToolIOValue(section)], Path(output[0]))

rule typing_dump_summary_info:
    """
    Dumps the summary information in tabular format.
    """
    input:
        HITS_NUCL = rules.typing_add_allele_page_url.output.HITS_NUCL,
        HITS_PEPT = rules.typing_add_allele_page_url.output.HITS_PEPT,
        INFORMS_ST = rules.typing_detect_sequence_type.output.INFORMS
    output:
        FILE = 'typing/{scheme}/summary_out.{ext}' # sequence_typing.OUTPUT_SUMMARY
    params:
        scheme_name = lambda wildcards: wildcards.scheme,
        ext = lambda wildcards: wildcards.ext,
        include_hashing = lambda wildcards: config['sequence_typing'][wildcards.scheme].get('include_hashing')
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
        hits = snakemakeutils.load_object(Path(input.HITS_NUCL)) + snakemakeutils.load_object(Path(input.HITS_PEPT))
        for io in hits:
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
