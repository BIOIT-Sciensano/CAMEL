from pathlib import Path

import pandas as pd

from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camelcore.app.io.tooliofile import ToolIOFile
from camel.app.core.snakemake.step import Step
from camel.app.core.snakemake import snakemakeutils


rule straingst_kmerize:
    """
    Runs straingst kmerize on reads.
    """
    input:
        FASTQ = 'fq_dict.io'
    output:
        HDF5 = 'straingst/{read_type}/tool/hdf5-straingst.io'
    params:
        dir_ = lambda wildcards: f'straingst/{wildcards.read_type}/tool',
        read_type = lambda wildcards: wildcards.read_type
    run:
        from camel.app.tools.strainge.straingstkmerize import StrainGSTKmerize
        straingst_kmerize_ = StrainGSTKmerize()
        fq_in = FastqInput.from_fq_dict(Path(input.FASTQ), f'{params.read_type}')
        if params.read_type == 'ont':
            straingst_kmerize_.add_input_files({'FASTQ': [fq_in.se[0]]})
        if params.read_type == 'illumina':
            straingst_kmerize_.add_input_files({'FASTQ': [fq_in.pe[0], fq_in.pe[1]]})
        step = Step(str(rule), straingst_kmerize_, dir_=Path(str(params.dir_)))
        step.run_step()
        snakemakeutils.dump_io_outputs(straingst_kmerize_, output)

rule straingst_run:
    """
    Runs straingst run on the generated HDF5 against the selected database.
    """
    input:
        HDF5 = rules.straingst_kmerize.output.HDF5
    output:
        TSV_STATS = 'straingst/{read_type}/tool/straingst-stats.io',
        TSV_STRAINS = 'straingst/{read_type}/tool/straingst-strains.io',
        INFORMS = 'straingst/{read_type}/tool/informs.io'
    params:
        dir_ = lambda wildcards: f'straingst/{wildcards.read_type}/tool',
        straingst_db = config['straingst']['db']
    run:
        from camel.app.tools.strainge.straingstrun import StrainGSTRun
        straingst_run_ = StrainGSTRun()
        snakemakeutils.add_io_inputs(straingst_run_, input)
        straingst_run_.add_input_files({'DB_HDF5': [ToolIOFile(Path(params.straingst_db))]})
        step = Step(str(rule), straingst_run_, dir_=Path(str(params.dir_)))
        step.run_step()
        snakemakeutils.dump_io_outputs(straingst_run_, output)

rule straingst_report:
    """
    Creates the report section for strainGST.
    """
    input:
        TSV = rules.straingst_run.output.TSV_STRAINS,
        INFORMS_straingst = rules.straingst_run.output.INFORMS
    output:
        VAL_HTML = 'straingst/{read_type}/report/html.iob'
    params:
        dir_ = lambda wildcards: f'straingst/{wildcards.read_type}/report/',
        sample_name = config['input']['sample_name'],
        read_type = lambda wildcards: wildcards.read_type
    run:
        from camel.app.tools.strainge.straingstreporter import StrainGSTReporter
        straingst_reporter = StrainGSTReporter()
        snakemakeutils.add_io_inputs(straingst_reporter, input)
        straingst_reporter.update_parameters(sample_name=params.sample_name, suffix=str(params.read_type))
        step = Step(str(rule), straingst_reporter, dir_=Path(str(params.dir_)))
        step.run_step()
        snakemakeutils.dump_io_outputs(straingst_reporter, output)

rule straingst_report_empty:
    """
    Creates an empty HTML report for the StrainGST analysis.
    """
    output:
        VAL_HTML = 'straingst/{read_type}/report/html-empty.iob'
    params:
        dir_ = lambda wildcards: f'straingst/{wildcards.read_type}/report/'
    run:
        from camel.app.core.snakemake import snakepipelineutils
        from camel.app.tools.strainge.straingstreporter import StrainGSTReporter
        snakepipelineutils.create_empty_report_section(StrainGSTReporter.TITLE, Path(output.VAL_HTML))

rule straingst_dump_summary_info:
    """
    Dumps the summary information for the StrainGST workflow in tabular format.
    """
    input:
        TSV = rules.straingst_run.output.TSV_STRAINS
    output:
        FILE = 'straingst/{read_type}/summary/summary_out.{ext}'
    params:
        read_type = lambda wildcards: wildcards.read_type,
        ext = lambda wildcards: wildcards.ext
    run:
        tsv_straingst = snakemakeutils.load_object(Path(input.TSV))[0].path
        straingst_table = pd.read_table(tsv_straingst)
        data_out = [
            (f'straingst_closest_strain_{params.read_type}', str(straingst_table.iloc[0]["strain"])),
            (f'straingst_breadth_of_coverage_{params.read_type}', str(straingst_table.iloc[0]["cov"])),
            (f'straingst_evenness_{params.read_type}', str(straingst_table.iloc[0]["even"])),
            (f'straingst_relative_abundance_{params.read_type}', str(straingst_table.iloc[0]["rapct"])),
            (f'straingst_score_{params.read_type}', str(straingst_table.iloc[0]["score"]))
        ]
        snakemakeutils.export_summary(data_out, Path(output.FILE), str(params.ext), 'straingst')
