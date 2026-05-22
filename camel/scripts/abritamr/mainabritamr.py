#!/usr/bin/env python
import dataclasses
from pathlib import Path
from typing import Any

import click
import yaml
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.cli import cliutils
from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.model import BaseOptions, BaseOutput
from camel.scripts.abritamr import CONFIG_DATA, SNAKEFILE_MAIN
from camel.snakefiles import assembly

SUPPORTED_SPECIES = [
    'Acinetobacter_baumannii',
    'Burkholderia_cepacia',
    'Burkholderia_pseudomallei',
    'Burkholderia_mallei',
    'Campylobacter',
    'Citrobacter_freundii',
    'Clostridioides_difficile',
    'Corynebacterium_diphtheriae',
    'Enterobacter_asburiae',
    'Enterobacter_cloacae',
    'Enterococcus_faecalis',
    'Enterococcus_faecium',
    'Escherichia',
    'Klebsiella_oxytoca',
    'Klebsiella_pneumoniae',
    'Neisseria_gonorrhoeae',
    'Neisseria_meningitidis',
    'Pseudomonas_aeruginosa',
    'Salmonella',
    'Serratia_marcescens',
    'Staphylococcus_aureus',
    'Staphylococcus_pseudintermedius',
    'Streptococcus_agalactiae',
]


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for abritAMR.
    """
    species: str = dataclasses.field(metadata={'choices': SUPPORTED_SPECIES})
    threads: int = dataclasses.field(default=1, metadata={
        'help': 'Number of threads',
        'show_default': True})
    working_dir: Path = dataclasses.field(default=Path('working'), metadata={'help': 'Working directory'})


@dataclasses.dataclass(frozen=True)
class Output(BaseOutput):
    """
    Tool output.
    """
    output_html: Path = dataclasses.field(metadata={'help': 'Output HTML file'})
    output_dir: Path | None = dataclasses.field(metadata={'help': 'Output directory'})
    output_tsv: Path | None = dataclasses.field(default=None, metadata={'help': 'Output TSV file'})


class MainAbriTAMR(BaseScript[FastaInput, Output, Options]):
    """
    Main class to run the AbriTAMR standalone pipeline.
    """

    def __init__(
            self,
            in_: FastaInput,
            out: BaseOutput,
            opts: Options,
    ) -> None:
        """
        Initializes the main class.
        :param in_: Script input
        :param out: Script output
        :param opts: Script options
        :return: None
        """
        super().__init__(
            name='AbritAMR',
            title='AbritAMR',
            version='0.3.0',
            script_in=in_,
            script_out=out,
            script_opts=opts
        )

    def _execute(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        script_in: FastaInput = self._script_in.create_symlinks(self._script_opts.working_dir)
        config_file = self.__construct_config_file(script_in)
        snakepipelineutils.run_snakemake(
            snakefile=SNAKEFILE_MAIN,
            targets=[],
            working_dir=self._script_opts.working_dir,
            config_path=config_file,
            threads=self._script_opts.threads
        )

    def __construct_config_file(self, script_in: FastaInput) -> str:
        """
        Constructs the configuration file.
        :param script_in: Script input
        :return: Configuration file
        """
        config_data: dict[str, Any] = {
            'input': {
                'fasta': {
                    'name': script_in.name,
                    'path': str(script_in.fasta)},
                'type': model.InputType.FASTA.value,
                'sample_name': self._script_in.name,
            },
            'output': {
                'html': str(self._script_out.output_html.absolute()),
                'dir': str(self._script_out.output_dir.absolute()),
                'tsv': str(self._script_out.output_tsv.absolute()),
            },
            'working_dir': str(self._script_opts.working_dir.absolute()),
            'script_info': self.info()
        }
        path_assembly = self._script_opts.working_dir / assembly.OUTPUT_FASTA
        path_assembly.parent.mkdir(parents=True, exist_ok=True)
        snakemakeutils.dump_object([ToolIOFile(self._script_in.fasta)], path_assembly)
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        config_data['abritamr']['species'] = self._script_opts.species
        config_data['analyses_selected'] =['abritamr']
        return snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)


@click.command(name='abritamr', short_help='Wrapper for AbritAMR')
@cliutils.add_click_options_from_dataclass(FastaInput)
@cliutils.add_click_options_from_dataclass(Output)
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for AbritAMR.
    """
    script = MainAbriTAMR(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=Output(**cliutils.from_kwargs(Output, kwargs)),
        opts=Options(**cliutils.from_kwargs(Options, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
