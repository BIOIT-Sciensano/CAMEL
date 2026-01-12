import dataclasses
import json
from pathlib import Path
from typing import Any

import click

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.cameltesthelper import running_dir
from camel.app.core.reports import reportutils
from camel.app.core.snakemake.snakemakeutils import IOEncoder, io_hook
from camel.app.core.utils import fastqutils, fastautils
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.inputhelper import helper_by_input_type
from camel.app.scriptutils.inputhelper.illuminahelper import (
    IlluminaTrimmingOpts,
    IlluminaAssemblyOpts,
)
from camel.app.scriptutils.inputhelper.inputhelperbase import TrimmingOpts, AssemblyOpts
from camel.app.scriptutils.inputhelper.onthelper import ONTTrimmingOpts, ONTAssemblyOpts


@click.command()
@basescriptutils.add_input_opts()
@click.option('--working-dir', type=click.Path(exists=True, path_type=Path), default=Path.cwd(), help='Working directory')
@click.option('--prepare', type=click.Choice(['fasta', 'fastq']), required=True, help='Target file type')
@click.option('--out', type=click.Path(path_type=Path), required=False, help='Output file')
@click.option('--threads', type=int, default=1, help='Nb. of threads')
@cliutils.add_click_options_from_dataclass(IlluminaTrimmingOpts)
@cliutils.add_click_options_from_dataclass(ONTTrimmingOpts, skip=[f.name for f in dataclasses.fields(TrimmingOpts)])
@cliutils.add_click_options_from_dataclass(IlluminaAssemblyOpts)
@cliutils.add_click_options_from_dataclass(ONTAssemblyOpts, skip=[f.name for f in dataclasses.fields(AssemblyOpts)])
def cli_test_helper(**kwargs: dict[str, Any]) -> None:
    """
    Test CLI to test the helper and argument parsing.
    :param kwargs: Command line arguments
    :return: None
    """
    # Setup report and script input
    report = reportutils.init_report(Path(str(kwargs['working_dir']), 'report.html'), key='Test report', title='Test report')
    script_input = basescriptutils.parse_script_input(kwargs)

    # Initialize the helper and prepare the input
    helper = helper_by_input_type[script_input.type_](
        dir_=kwargs['working_dir'], name=script_input.name
    )
    helper.set_opts(*helper.opts_from_cli(kwargs))

    # Collect the output
    data_out = {}
    if kwargs['prepare'] == 'fasta':
        fasta = helper.prepare_fasta_input(script_input, report)
        data_out['fasta'] = str(fasta)
    else:
        fastq = helper.prepare_fastq_input(script_input, report)
        data_out['fastq'] = fastq.to_fq_dict()
    data_out['informs'] = helper.informs

    with Path(str(kwargs['out'])).open('w') as handle:
        json.dump(data_out, handle, indent=2, cls=IOEncoder)

base_opts_ont = [
    '--fastq-se', str(Path(config.dir_testdata, 'input_helper/reads_ont.fastq.gz')),
    '--input-type', 'ont',
    '--working-dir', str(running_dir),
    '--threads', '4'
]
base_opts_illumina = [
    '--fastq-pe',
    str(Path(config.dir_testdata, 'input_helper/reads_ilmn_1.fastq.gz')),
    str(Path(config.dir_testdata, 'input_helper/reads_ilmn_2.fastq.gz')),
    '--input-type', 'illumina',
    '--threads', '4',
]

def test_ont_helper_fastq_trim(running_dir: Path) -> None:
    """
    Tests the ONT helper with read trimming.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    min_qual = 12
    min_len = 750
    result = cliutils.invoke(
        cli_test_helper,[
            *base_opts_ont,
            '--working-dir', str(running_dir),
            '--prepare', 'fastq',
            '--trim-reads',
            '--ont-min-len', str(min_len),
            '--ont-min-qual', str(min_qual),
            '--out', str(path_out)
        ],
    )
    # Check the output
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the reads were trimmed
    bases_in = fastqutils.count_bases(Path(config.dir_testdata, 'input_helper/reads_ont.fastq.gz'))
    bases_out = fastqutils.count_bases(data_out['fastq']['SE'][0].path)
    assert bases_in > bases_out, "Reads were not trimmed"

    # Checks if the correct parameters were used
    assert data_out['informs'][0]['min_length'] == min_len, "Wrong minimum read length"
    assert data_out['informs'][0]['min_qual'] == min_qual, "Wrong minimum quality"

def test_ont_helper_fastq_no_trim(running_dir: Path) -> None:
    """
    Tests the ONT helper with no read trimming.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_ont,
        '--working-dir', str(running_dir),
        '--prepare', 'fastq',
        '--out', path_out,
    ])
    # Check the output
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check that the reads were not trimmed
    bases_in = fastqutils.count_bases(Path(config.dir_testdata, 'input_helper/reads_ont.fastq.gz'))
    bases_out = fastqutils.count_bases(data_out['fastq']['SE'][0].path)
    assert bases_in == bases_out, "Reads were trimmed"

def test_illumina_helper_fastq_trim(running_dir: Path) -> None:
    """
    Tests the Illumina helper with read trimming.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_illumina,
        '--working-dir', str(running_dir),
        '--prepare', 'fastq',
        '--trim-reads',
        '--out', path_out,
    ])
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the reads were trimmed
    for idx, ori in enumerate(['1', '2']):
        bases_in = fastqutils.count_bases(Path(config.dir_testdata, f'input_helper/reads_ilmn_{ori}.fastq.gz'))
        bases_out = fastqutils.count_bases(data_out['fastq']['PE'][idx].path)
        assert bases_in > bases_out, "Reads were not trimmed"

    # Checks if the correct parameters were used
    assert len(data_out['informs']) == 1, "No trimming informs found"

def test_illumina_helper_fastq_no_trim(running_dir: Path) -> None:
    """
    Tests the Illumina helper without read trimming.
    :param running_dir: Running directory
    :return: None
    """
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_illumina,
        '--working-dir', str(running_dir),
        '--prepare', 'fastq',
    ])
    assert result.exit_code == 0

def test_ont_helper_fasta_trim(running_dir: Path) -> None:
    """
    Tests the ONT helper with read trimming and FASTA output.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    min_qual = 12
    min_len = 750
    result = cliutils.invoke(
        cli_test_helper,[
            *base_opts_ont,
            '--working-dir', str(running_dir),
            '--prepare', 'fasta',
            '--trim-reads',
            '--ont-min-len', str(min_len),
            '--ont-min-qual', str(min_qual),
            '--out', path_out,
        ],
    )
    # Check the output
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the assembly succeeded
    assert fastautils.count_reads(Path(data_out['fasta'])), "Assembly failed"

def test_ont_helper_fasta_no_trim(running_dir: Path) -> None:
    """
    Tests the ONT helper with no read trimming and FASTA output.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_ont,
        '--working-dir', str(running_dir),
        '--prepare', 'fasta',
        '--out', path_out,
    ])
    # Check the output
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the assembly succeeded
    assert fastautils.count_reads(Path(data_out['fasta'])), "Assembly failed"

def test_illumina_helper_fasta_trim(running_dir: Path) -> None:
    """
    Tests the Illumina helper with read trimming and FASTA output.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_illumina,
        '--working-dir', str(running_dir),
        '--prepare', 'fasta',
        '--trim-reads',
        '--out', path_out,
    ])
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the assembly succeeded
    assert fastautils.count_reads(Path(data_out['fasta'])), "Assembly failed"

def test_illumina_helper_fasta_no_trim(running_dir: Path) -> None:
    """
    Tests the Illumina helper without read trimming and FASTA output.
    :param running_dir: Running directory
    :return: None
    """
    path_out = running_dir / 'helper_out.json'
    result = cliutils.invoke(cli_test_helper, [
        *base_opts_illumina,
        '--working-dir', str(running_dir),
        '--prepare', 'fasta',
        '--out', path_out,
    ])
    assert result.exit_code == 0
    with path_out.open() as handle:
        data_out = json.load(handle, object_hook=io_hook)

    # Check if the assembly succeeded
    assert fastautils.count_reads(Path(data_out['fasta'])), "Assembly failed"
