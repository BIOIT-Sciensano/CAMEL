#!/usr/bin/env python
import dataclasses
from importlib.resources import files
from typing import Any

import click
import yaml

from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.mockpipeline import SNAKEFILE_MAIN

CUSTOM_ANALYSES = [
    'human_read_scrubbing',
    'kraken2',
    'confindr',
    'ncbi_amr',
    'snpit',
]


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """
    analyses: list[str] = dataclasses.field(default_factory=list)


class MainMockPipeline(BasePipe):
    """
    Base-class for the mock pipeline.
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
        super().__init__(
            name='Mock pipeline',
            version='1.0',
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
        # Parse template data
        with open(str(files('camel').joinpath('scripts/mockpipeline/config_data.yml'))) as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            DB_ROOT=config.dir_db
        )
        data_template: dict[str, Any] = yaml.safe_load(yaml_text)
        self.check_dbs(data_template)

        self._script_out.dir.mkdir(parents=True, exist_ok=True)

        # Add the base config data
        config_data = self.get_config_data()
        config_data['analyses'] = self._opts_custom.analyses
        basepipeutils.dict_merge(config_data, data_template)
        path_config = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run the Snakefile
        snakepipelineutils.run_snakemake(
            snakefile=self._snakefile,
            config_path=path_config,
            targets=[],
            working_dir=self._script_opts.working_dir,
            threads=self._script_opts.threads)
        self._export_assembly()


@click.command(name='mock_pipeline', short_help='Dummy pipeline for testing')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option('--analyses', type=str, help=f"Comma-separated list of analyses to run ({', '.join(CUSTOM_ANALYSES)})")
def main(**kwargs) -> None:
    """
    Dummy pipeline for testing.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
    )
    mock_pipe = MainMockPipeline(script_input, script_out, script_opts, custom_opts)
    mock_pipe.run()


if __name__ == '__main__':
    initialize_logging()
    main()
