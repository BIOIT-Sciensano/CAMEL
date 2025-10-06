from pathlib import Path

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.pipeline.step import Step
from camel.app.snakemake import snakemakeutils
from camel.resources.snakefile import assembly


rule quast_quast:
    """
    Runs quast on the assembly.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA,
        IO = 'fq_dict.io'
    output:
        TSV = 'quast/output/tsv.io',
        HTML = 'quast/output/html.iob', # quast.OUTPUT_REPORT
        DIR = 'quast/output/dir.io',
        INFORMS = 'quast/output/informs.io' # quast.OUTPUT_INFORMS
    params:
        dir_ = 'quast/output',
        input_type = config['input_type'],
        fasta = config.get('reference', {}).get('fasta'),
        gff = config.get('reference', {}).get('gff3')
    run:
        from camel.app.io.tooliofile import ToolIOFile
        from camel.app.tools.quast.quast import Quast

        # Create output directory
        dir_out = Path(params.dir_) / 'quast_out'
        dir_out.mkdir(exist_ok=True, parents=True)

        # Create tool
        quast_ = Quast()

        # Add input
        snakemakeutils.add_pickle_inputs(quast_, input, excluded_keys=['IO'])
        fq_dict = snakemakeutils.load_object(Path(input.IO))
        if params.input_type in ('illumina', 'hybrid'):
            quast_.add_input_files({'FASTQ_PE': fq_dict['PE']})
        if params.input_type in ('ont', 'hybrid'):
            quast_.add_input_files({'FASTQ_nanopore': fq_dict['SE']})

        # Add reference genome files (if available)
        if (params.fasta is not None) and (params.gff is not None):
            quast_.add_input_files({
                'FASTA_Ref': [ToolIOFile(Path(params.fasta))],
                'GFF3_Ref': [ToolIOFile(Path(params.gff))],
            })
        else:
            logger.warning(f'No reference genome provided, skipping analysis for QUAST')

        # Run tool
        quast_.update_parameters(conserved_genes_finding=False)
        step = Step(rule_name=str(rule), tool=quast_, dir_=dir_out)
        step.run()

        # Collect output
        snakemakeutils.dump_tool_outputs(quast_, output, ignore_missing_output=True)
        snakemakeutils.dump_object([ToolIODirectory(dir_out)], Path(output.DIR))

rule quast_busco:
    """
    Runs BUSCO on the assembly to check completeness.
    BUSCO is ran outside of QUAST because of dependency issues with the BUSCO installation bundled with QUAST.
    """
    input:
        FASTA = assembly.OUTPUT_FASTA
    output:
        TXT = 'quast/busco/txt.io',
        INFORMS = 'quast/busco/informs.io' # quast.OUTPUT_INFORMS_BUSCO
    params:
        dir_ = 'quast/busco',
        lineage_dataset = 'bacteria_odb10'
    threads: 8
    run:
        from camel.app.tools.busco.busco import Busco
        busco = Busco()
        snakemakeutils.add_pickle_inputs(busco, input)
        step = Step(rule_name=str(rule), tool=busco, dir_=Path(params.dir_))
        busco.update_parameters(lineage_dataset=params.lineage_dataset, threads=str(threads))
        step.run()
        snakemakeutils.dump_tool_outputs(busco, output)

rule quast_report:
    """
    Creates a report for the quast workflow.
    """
    input:
        TSV = rules.quast_quast.output.TSV,
        HTML = rules.quast_quast.output.HTML,
        FASTA =  assembly.OUTPUT_FASTA,
        DIR = rules.quast_quast.output.DIR,
        INFORMS_quast = rules.quast_quast.output.INFORMS,
        INFORMS_assembler = assembly.get_command_informs(config),
        INFORMS_busco = rules.quast_busco.output.INFORMS
    output:
        HTML = 'quast/report/html.iob' # quast.OUTPUT_REPORT
    params:
        dir_ = 'quast/report',
        name = config['sample_name']
    run:
        from camel.app.tools.quast.quastreporter import QuastReporter
        reporter = QuastReporter()
        reporter.update_parameters(name=params.name)
        snakemakeutils.add_pickle_inputs(reporter, input, excluded_keys=['INFORMS_assembler'])
        reporter.add_input_informs({
            'assembler': ', '.join(snakemakeutils.load_object(Path(x))['_name'] for x in input.INFORMS_assembler)
            if input.INFORMS_assembler else 'n/a'
        })
        step = Step(rule_name=str(rule), tool=reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(reporter, output)

rule quast_create_summary_out:
    """
    Creates the tabular summary output for QUAST.
    """
    input:
        INFORMS = assembly.get_command_informs(config),
        INFORMS_filtering = 'assembly/filtering/informs.io',
        TSV = rules.quast_quast.output.TSV
    output:
        FILE = 'quast/summary/summary_quast.{ext}' # quast.OUTPUT_SUMMARY
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        keys_kept = [
            {'key': '# contigs', 'name': 'nb_contigs'},
            {'key': 'Total length', 'name': 'total_length'},
            {'key': 'Reference length', 'name': 'total_length_ref'},
            {'key': 'N50', 'name': 'n50'},
            {'key': 'Genome fraction (%)', 'name': 'genome_fraction'},
            {'key': 'Duplication ratio', 'name': 'dupl_ratio'},
            {'key': 'Avg. coverage depth', 'name': 'avg_coverage'},
            {'key': 'Reference avg. coverage depth', 'name': 'avg_coverage_ref'},
            {'key': 'Coverage >= 1x (%)', 'name': 'positions_covered_1x'},
            {'key': 'Reference coverage >= 1x (%)', 'name': 'positions_covered_1x_ref'},
            {'key': 'tool_versions', 'name': 'tool_versions'},
            {'key': 'filtering_tool_version', 'name': 'filtering_tool_version'}
        ]

        # Parse QUAST report
        path_tsv = snakemakeutils.load_object(Path(input.TSV))[0].path
        data_quast = {}
        with path_tsv.open() as handle:
            for line in handle.readlines():
                key, value = line.strip().split('\t')
                data_quast[key] = value

        # Add assembler version
        tool_names = [snakemakeutils.load_object(Path(x))['_name'] for x in input.INFORMS] if input.INFORMS else ['n/a']
        data_quast['tool_version'] = ', '.join(tool_names)
        data_quast['filtering_tool_version'] = snakemakeutils.load_object(Path(input.INFORMS_filtering))['_name']

        # Create output
        data_out = [(f"assembly_{row['name']}", data_quast.get(row['key'], '-')) for row in keys_kept]
        snakemakeutils.export_summary(data_out, Path(output.FILE), str(params.ext), 'quast')
