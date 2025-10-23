from pathlib import Path

import pandas as pd

from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils
from camel.scripts.bacilluspipeline.snakefile import btyper as bt


rule btyper_run:
    """
    Runs btyper on assembled contigs.
    """
    input:
        FASTA = bt.INPUT_FASTA
    output:
        TSV = 'btyper/tool/val-btyper.io',
        INFORMS = 'btyper/tool/informs.io'
    params:
        dir_ = 'btyper/tool'
    run:
        from camel.app.tools.btyper.btyper import BTyper
        btyper = BTyper()
        snakemakeutils.add_pickle_inputs(btyper, input)
        step = Step(rule_name=str(rule), tool=btyper, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(btyper, output)

rule btyper_report:
    """
    Creates the report section for btyper.
    """
    input:
        TSV = rules.btyper_run.output.TSV,
        INFORMS_btyper = rules.btyper_run.output.INFORMS
    output:
        HTML = 'btyper/report/html.iob'
    params:
        dir_ = 'btyper/report',
        sample_name = config['sample_name']
    run:
        from camel.app.tools.btyper.btyperreporter import BTyperReporter
        btyper_reporter = BTyperReporter()
        snakemakeutils.add_pickle_inputs(btyper_reporter, input)
        btyper_reporter.update_parameters(sample_name=params.sample_name)
        step = Step(rule_name=str(rule), tool=btyper_reporter, dir_=Path(params.dir_))
        step.run()
        snakemakeutils.dump_tool_outputs(btyper_reporter, output)

rule btyper_report_empty:
    """
    Creates an empty HTML report for the BTyper analysis.
    """
    output:
        HTML = 'btyper/report/html-empty.iob'
    params:
        dir_ = 'btyper/report'
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.btyper.btyperreporter import BTyperReporter
        snakepipelineutils.create_empty_report_section(BTyperReporter.TITLE, Path(output.HTML))

rule btyper_dump_summary_info:
    """
    Dumps the summary information for the BTyper workflow in tabular format.
    """
    input:
        TSV = rules.btyper_run.output.TSV
    output:
        FILE = 'btyper/summary/summary_out.{ext}'
    params:
        ext = lambda wildcards: wildcards.ext
    run:
        tsv_btyper = snakemakeutils.load_object(Path(input.TSV))[0].path
        btyper_table = pd.read_table(tsv_btyper)

        # Parse virulence and capsule genes
        subtable_virulence = btyper_table.iloc[:, 6:12]
        virulence_data = [f"{col}:{subtable_virulence.iloc[0][col]}" for col in subtable_virulence.columns]
        subtable_capsule = btyper_table.iloc[:,12:15]
        capsule_data = [f"{col}:{subtable_capsule.iloc[0][col]}" for col in subtable_capsule.columns]

        # Create the summary data
        data_summary = [
            ('btyper_closest_species_ANI', btyper_table['species(ANI)'][0]),
            ('btyper_panC_type', btyper_table['Adjusted_panC_Group(predicted_species)'][0]),
            ('btyper_virulence_genes', ';'.join(virulence_data)),
            ('btyper_capsule_genes', ';'.join(capsule_data)),
            ('btyper_bt_genes', btyper_table['Bt(genes)'][0]),
        ]
        snakemakeutils.export_summary(data_summary, Path(output.FILE), str(params.ext), 'btyper')
