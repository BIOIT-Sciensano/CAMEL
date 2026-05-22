from pathlib import Path

from camelcore.app.io.tooliofile import ToolIOFile
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.bacilluspipeline.snakefile import ani


rule fastani_run:
    """
    Runs fastani on assembled contigs.
    """
    input:
        FASTA_Q = ani.INPUT_FASTA
    output:
        TSV = 'ani/tool/val-ani.io',
        INFORMS = 'ani/tool/informs.io'
    params:
        fastani_ref = config['fastani']['path'],
        dir_ = 'ani/tool'
    run:
        from camel.app.tools.fastani.fastani import FastANI
        fastani = FastANI()
        snakemakeutils.add_io_inputs(fastani, input)
        fastani.add_input_files({'TSV_FASTA_R': [ToolIOFile(Path(params.fastani_ref))]})
        step = Step(rule_name=str(rule), tool=fastani, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_io_outputs(fastani, output)

rule fastani_report:
    """
    Creates the report section for fastani.
    """
    input:
        TSV = rules.fastani_run.output.TSV,
        INFORMS_fastani = rules.fastani_run.output.INFORMS
    output:
        HTML = 'ani/report/html.iob'
    params:
        dir_ = 'ani/report',
        sample_name = config['input']['sample_name'],
        species = config['species']
    run:
        from camel.app.tools.fastani.fastanireporter import FastANIReporter
        ani_report = FastANIReporter()
        snakemakeutils.add_io_inputs(ani_report, input)
        ani_report.update_parameters(sample_name=params.sample_name, species=params.species)
        step = Step(rule_name=str(rule), tool=ani_report, dir_=Path(str(params.dir_)))
        step.run()
        snakemakeutils.dump_io_outputs(ani_report, output)

rule fastani_report_empty:
    """
    Creates an empty HTML report for the FastANI analysis.
    """
    output:
        HTML = 'ani/report/html-empty.iob'
    params:
        dir_ = 'ani/report'
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.fastani.fastanireporter import FastANIReporter
        snakepipelineutils.create_empty_report_section(FastANIReporter.TITLE, Path(output.HTML))

rule fastani_dump_summary_info:
    """
    Dumps the summary information for the FastANI workflow in tabular format.
    """
    input:
        TSV = rules.fastani_run.output.TSV
    output:
        TSV = 'ani/summary/summary_out.{ext}'
    run:
        import pandas as pd
        tsv_fastani = snakemakeutils.load_object(Path(input.TSV))[0].path
        fastani_table = pd.read_table(tsv_fastani, header=None)
        with open(output.TSV, 'w') as handle:
            handle.write('fastani_closest_species\t{}'.format([Path(fastani_table[1][0]).stem, str(fastani_table[2][0])+'%']))
            handle.write('\n')
