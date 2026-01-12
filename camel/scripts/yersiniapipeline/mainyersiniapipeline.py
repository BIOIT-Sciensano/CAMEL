#!/usr/bin/env python
import dataclasses
from pathlib import Path

import click
import yaml

from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.scripts.yersiniapipeline import SNAKEFILE_MAIN, CONFIG_DATA

CUSTOM_ANALYSES = [
    'ampc',
    'amrfinder',
    'cgmlst',
    'cgmlst_enterobase',
    'cgmlst_ye',
    'cgmlst_yp',
    'confindr',
    'human_read_scrubbing',
    'kraken2',
    'mlst',
    'mlst_mcnally',
    'mob_suite',
    'resfinder4',
    'rmlst',
    'vfdb_core',
]


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """
    analyses: list[str] = dataclasses.field(default_factory=list)


class MainYersiniaPipeline(BasePipe):
    """
    Main class to run the Yersinia pipeline.
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
            name='Yersinia pipeline',
            title='<i>Yersinia</i> pipeline',
            version='1.2',
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
        with open(CONFIG_DATA) as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            DB_ROOT=config.dir_db,
            QC_SCHEME='cgmlst' if 'cgmlst' in  self._opts_custom.analyses else 'mlst',
            K2_DB=str({
                'small': Path(config.dir_db, '/scratch/bebog/camel_dbs/kraken2/k2_standard_08_20251015'), # TODO FIX!
                'full': Path(config.dir_db, 'kraken2_microbial', 'latest')
            }['small' if self._script_opts.kraken2_small_db else 'full'])
        )
        data_template = yaml.safe_load(yaml_text)
        self.check_dbs(data_template)
        self._script_out.dir.mkdir(parents=True, exist_ok=True)

        # Add the base config data
        config_data = self.get_config_data()
        config_data['analyses'] = self._opts_custom.analyses
        if 'cgmlst' in self._opts_custom.analyses:
            config_data['analyses'].append('species')
        config_data['sequence_typing'] = {'options': {'method': self._script_opts.typing_method}}
        config_data['gene_detection'] = {'options': {'method': self._script_opts.gene_detection_method}}
        basepipeutils.dict_merge(config_data, data_template)
        path_config = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run the Snakefile
        self.run_snakefile(path_config)
        self._export_assembly()


@click.command(name='yersinia_pipeline', short_help='Pipeline for the complete characterization of Yersinia isolates')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option('--analyses', type=str, help=f"Comma-separated list of analyses to run ({', '.join(CUSTOM_ANALYSES)})")
def main(**kwargs) -> None:
    """
    Pipeline for the complete characterization of Yersinia isolates.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
    )
    pipeline = MainYersiniaPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.prepare_input()
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
