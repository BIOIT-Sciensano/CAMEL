#!/usr/bin/env python
import dataclasses

import click
import yaml

from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.neisseriapipeline import SNAKEFILE_MAIN, CONFIG_DATA

CUSTOM_ANALYSES = [
    'amrfinder',
    'bast',
    'cgmlst',
    'cgmlst_v3',
    'confindr',
    'feta',
    'fhbp',
    'gmats',
    'human_read_scrubbing',
    'kraken2',
    'mendevar',
    'mlst',
    'pora',
    'porb',
    'resfinder4',
    'resistance_genes',
    'rmlst',
    'rplf',
    'serogroup',
    'vaccine_targets',
    'variant_calling'
]


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
            opts_custom: Options
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
            snakefile=SNAKEFILE_MAIN
        )
        self._opts_custom = opts_custom

    def _execute(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        # Parse template data
        with CONFIG_DATA.open() as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            QC_SCHEME='cgmlst' if 'cgmlst' in self._opts_custom.analyses else 'mlst',
        )
        data_template = yaml.safe_load(yaml_text)

        # Add the base config data
        config_data = self.get_config_data()
        basepipeutils.dict_merge(config_data, data_template)
        config_data['analyses'] = self._opts_custom.analyses
        config_data['sequence_typing']['options'] = {'method': self._script_opts.typing_method}
        path_config = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run the Snakefile
        self.run_snakefile(path_config)


@click.command(name='neisseria_pipeline', short_help='Pipeline for the complete characterization of Neisseria meningitidis isolates')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option(
    '--analyses',
    type=str,
    help=f"Comma-separated list of analyses to run: {', '.join(CUSTOM_ANALYSES)}")
def main(**kwargs) -> None:
    """
    Pipeline for the complete characterization of Neisseria meningitidis isolates.
    """
    BasePipe.check_analyses_option(kwargs['analyses'], CUSTOM_ANALYSES)
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [])
    pipeline = MainNeisseriaPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.prepare_input()
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
