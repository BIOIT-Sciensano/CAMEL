#!/usr/bin/env python
import dataclasses

import click
import yaml

from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.neisseriapipeline import SNAKEFILE_MAIN, CONFIG_DATA


@dataclasses.dataclass(frozen=True)
class Options:
    """
    Pipeline-specific options.
    """

    analyses: list[str] = dataclasses.field(default_factory=list)


class MainNeisseriaPipeline(BasePipe):
    """
    Main class to run the Neisseria pipeline.
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
        """
        super().__init__(
            name='Neisseria pipeline',
            title='<i>Neisseria</i> pipeline',
            version='1.5.0',
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
        with CONFIG_DATA.open() as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            DB_ROOT=config.dir_db,
            QC_SCHEME='cgmlst' if 'cgmlst' in self._opts_custom.analyses else 'mlst',
        )
        data_template = yaml.safe_load(yaml_text)

        # Add the base config data
        config_data = self.get_config_data()
        basepipeutils.dict_merge(config_data, data_template)
        config_data['analyses_selected'] = self._opts_custom.analyses
        config_data['sequence_typing']['options'] = {
            'method': self._script_opts.typing_method
        }
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

        # Create the config file and run snakefile
        self._script_out.dir.mkdir(parents=True, exist_ok=True)
        path_config = snakepipelineutils.generate_config_file(
            config_data, self._script_opts.working_dir
        )
        self.run_snakefile(path_config)

        # Additional export for the assembly
        self._export_assembly()


@click.command(
    name='neisseria_pipeline',
    short_help='Pipeline for the complete characterization of Neisseria meningitidis isolates',
)
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option(
    '--analyses',
    type=str,
    help=f"Comma-separated list of analyses to run ({', '.join(basepipeutils.get_custom_analyses(CONFIG_DATA))})",
)
def main(**kwargs) -> None:
    """
    Pipeline for the complete characterization of Neisseria meningitidis isolates.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_analyses = basepipeutils.get_custom_analyses(CONFIG_DATA)
    BasePipe.check_analyses_option(kwargs['analyses'], custom_analyses)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else []
    )
    pipeline = MainNeisseriaPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.prepare_input()
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
