#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path

import click
import yaml

from camel.app.cli import cliutils
from camel.app.core.snakemake import snakemakeutils, snakepipelineutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import (
    NcbiHumanReadScrubber,
)
from camel.scripts.ncbihumanreadscrubber import CONFIG_DATA, SNAKEFILE_MAIN
from camel.snakefiles import human_read_scrubbing
from camel.version import __VERSION__


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the NCBI human read scrubber tool.
    """
    export_removed_reads: bool = dataclasses.field(default=False, metadata={
        'help': 'Export the removed reads'})


class MainNcbiHumanReadScrubber(BasePipe):
    """
    Main class to run the NCBI human read scrubber tool.
    """

    def __init__(
        self,
        in_: ScriptInput,
        out: ScriptOutput,
        opts: ScriptOptions,
        opts_custom: Options
    ) -> None:
        """
        Initializes the main class.
        :param in_: Script input
        :param out: Script output
        :param opts: General pipeline options
        :param opts_custom: Pipeline-specific options
        :return: None
        """
        tool_version = NcbiHumanReadScrubber().version
        super().__init__(
            name='NCBI human read scrubbing',
            title='NCBI human read scrubbing',
            version=f'{tool_version}+CAMEL_{__VERSION__}',
            script_in=in_,
            script_out=out,
            opts=opts,
            snakefile=SNAKEFILE_MAIN
        )
        self._opts_custom = opts_custom

    def _execute(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        config_file = self.__construct_config_file()
        self.run_snakefile(config_file)
        self._copy_output_files()

    def __construct_config_file(self) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_config_data()
        config_data['analyses_selected'] =['human_read_scrubbing']
        config_data['read_scrubbing'] = {}
        if self._opts_custom.export_removed_reads:
            config_data['read_scrubbing']['export_removed_reads'] = True
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        self._config_data = config_data
        return snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

    def _copy_output_files(self) -> None:
        """
        Copies the output files to the output directory.
        :return: None
        """
        # Copy the scrubbed reads
        output_files = snakemakeutils.load_object(self._script_opts.working_dir / human_read_scrubbing.get_output_io(
            self._config_data))
        if self._script_in.type_ is model.InputType.ILLUMINA:
            shutil.copyfile(output_files[0].path, Path(self._script_out.dir / f'{self._script_in.name}-scrubbed_R1.fastq.gz'))
            shutil.copyfile(output_files[1].path, Path(self._script_out.dir / f'{self._script_in.name}-scrubbed_R2.fastq.gz'))
        elif self._script_in.type_ is model.InputType.ONT:
            shutil.copyfile(output_files[0].path, Path(self._script_out.dir / f'{self._script_in.name}-scrubbed.fastq.gz'))
        elif self._script_in.type_ is model.InputType.FASTA:
            shutil.copyfile(output_files[0].path, Path(self._script_out.dir / f'{self._script_in.name}-scrubbed.fasta'))
        else:
            raise ValueError(f"Invalid input type: {self._script_in.type_}")

        # Copy the removed reads
        if not self._opts_custom.export_removed_reads:
            return
        if self._script_in.type_ is model.InputType.ONT:
            fq_removed = snakemakeutils.load_object(self._script_opts.working_dir / human_read_scrubbing.get_removed('fastq_se'))
            if len(fq_removed) == 0:
                logger.warning('No removed reads found')
                return
            shutil.copyfile(fq_removed[0].path, self._script_out.dir / f'{self._script_in.name}-removed.fastq.gz')
        elif self._script_in.type_ is model.InputType.ILLUMINA:
            fq_removed = snakemakeutils.load_object(self._script_opts.working_dir / human_read_scrubbing.get_removed('fastq_pe'))
            if len(fq_removed) == 0:
                logger.warning('No removed reads found')
                return
            shutil.copyfile(fq_removed[0].path, self._script_out.dir / f'{self._script_in.name}-removed_R1.fastq.gz')
            shutil.copyfile(fq_removed[1].path, self._script_out.dir / f'{self._script_in.name}-removed_R2.fastq.gz')
        elif self._script_in.type_ is model.InputType.FASTA:
            fa_removed = snakemakeutils.load_object(
                self._script_opts.working_dir / human_read_scrubbing.get_removed('fasta'))
            if len(fa_removed) == 0:
                logger.warning('No removed reads found')
                return
            shutil.copyfile(fa_removed[0].path, self._script_out.dir / f'{self._script_in.name}-removed.fasta')

@click.command(name='ncbi_human_read_scrubber', short_help='Wrapper for NCBI human read scrubber')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Runs the main script.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(**cliutils.from_kwargs(Options, kwargs))
    pipe_script = MainNcbiHumanReadScrubber(script_input, script_out, script_opts, custom_opts)
    pipe_script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
