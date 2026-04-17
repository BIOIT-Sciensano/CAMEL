#!/usr/bin/env python
import dataclasses
from pathlib import Path

import click
import yaml

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.enterococcuspipeline import SNAKEFILE_MAIN, CONFIG_DATA


DISABLED_FOR_SPP = ['mlst', 'cgmlst', 'variant_calling']


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """

    species: str = dataclasses.field(
        metadata={'choices': ['faecalis', 'faecium', 'spp']}
    )
    analyses: list[str] = dataclasses.field(default_factory=list)


class MainEnterococcusPipeline(BasePipe):
    """
    Main class to run the Enterococcus pipeline.
    """

    def __init__(
        self,
        in_: ScriptInput,
        out: ScriptOutput,
        opts: ScriptOptions,
        opts_custom: Options,
    ) -> None:
        """
        Initializes the main class.
        :param in_: Script input
        :param out: Script output
        :param opts: General pipeline options
        :param opts_custom: Pipeline-specific options
        :return: None
        """
        super().__init__(
            name='Enterococcus pipeline',
            title='<i>Enterococcus</i> pipeline',
            version='1.3.0',
            script_in=in_,
            script_out=out,
            opts=opts,
            snakefile=SNAKEFILE_MAIN,
        )
        self._opts_custom = opts_custom

    def _build_config(self) -> dict:
        """
        Builds the configuration data for Snakemake.
        :return: Configuration data
        """
        # Parse template data
        with open(CONFIG_DATA) as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            DB_ROOT=config.dir_db,
            QC_SCHEME='cgmlst' if 'cgmlst' in self._opts_custom.analyses else 'mlst',
            K2_DB=str(
                {
                    'small': Path(
                        config.dir_db,
                        '/scratch/bebog/camel_dbs/kraken2/k2_standard_08_20251015',
                    ),
                    # TODO FIX!
                    'full': Path(config.dir_db, 'kraken2_microbial', 'latest'),
                }['small' if self._script_opts.kraken2_small_db else 'full']
            ),
        )
        # Fill in placeholders
        data_template = yaml.safe_load(yaml_text)
        config_data = self.get_config_data()

        # Species information
        config_data['species'] = self._opts_custom.species
        config_data['species_name'] = data_template['species'][
            self._opts_custom.species
        ]['full_name']
        config_data['is_generic'] = self._opts_custom.species == 'spp'

        # Typing and gene detection
        config_data['sequence_typing'] = {
            'options': {'method': self._script_opts.typing_method}
        }
        config_data['gene_detection'] = {
            'options': {'method': self._script_opts.gene_detection_method}
        }

        # Merge with the template
        basepipeutils.dict_merge(config_data, data_template)

        # Analyses
        config_data['analyses_selected'] = self._opts_custom.analyses
        # Disable incompatible assays for generic species
        if self._opts_custom.species == 'spp':
            config_data['analyses_selected'] = [
                a for a in config_data['analyses'] if a not in DISABLED_FOR_SPP
            ]
            logger.warning(
                f"Generic 'Enterococcus' selected as species, disabling assays: {', '.join(DISABLED_FOR_SPP)}"
            )

        # Resolve species specific values
        config_data = basepipeutils.resolve_config(
            config_data, self._opts_custom.species
        )
        return config_data

    def _validate_config_data(self, config_data: dict) -> bool:
        """
        Validates the config data.
        :param config_data: Config data
        :return: True if valid, False otherwise
        """
        self.check_dbs(config_data)
        return True

    def _execute(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        # Build and validate the config file
        config_data = self._build_config()
        self._validate_config_data(config_data)

        # Create the config file and run the snakefile
        self._script_out.dir.mkdir(parents=True, exist_ok=True)
        path_config = snakepipelineutils.generate_config_file(
            config_data, self._script_opts.working_dir
        )
        self.run_snakefile(path_config)

        # Additional export for the assembly
        self._export_assembly()


@click.command(
    name='enterococcus_pipeline',
    short_help='Pipeline for the complete characterization of Enterococcus isolates',
)
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option(
    '--analyses',
    type=str,
    help=f"Comma-separated list of analyses to run ({', '.join(basepipeutils.get_custom_analyses(CONFIG_DATA))})",
)
@cliutils.add_click_options_from_dataclass(Options, skip=['analyses'])
def main(**kwargs) -> None:
    """
    Runs the main script.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
        species=kwargs['species'],
    )
    pipe_script = MainEnterococcusPipeline(
        script_input, script_out, script_opts, custom_opts
    )
    pipe_script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
