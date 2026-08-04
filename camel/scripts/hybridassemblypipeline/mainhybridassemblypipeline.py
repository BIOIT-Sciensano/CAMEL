#!/usr/bin/env python
import dataclasses
from pathlib import Path

import click
import pandas as pd
import yaml

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.scripts.hybridassemblypipeline import CONFIG_DATA, SNAKEFILE_MAIN, TSV_BASECALLING_MODELS

BASECALLING_MODELS = [
    'r1041_e82_260bps_fast_g632',
    'r1041_e82_260bps_hac_g632',
    'r1041_e82_260bps_sup_g632',
    'r1041_e82_400bps_fast_g615',
    'r1041_e82_400bps_fast_g632',
    'r1041_e82_400bps_hac_g615',
    'r1041_e82_400bps_hac_g632',
    'r1041_e82_400bps_sup_g615',
    'r1041_e82_400bps_sup_v5.0.0',
    'r104_e81_hac_g5015',
    'r104_e81_sup_g5015',
    'r104_e81_sup_g610',
    'r941_min_hac_g507',
    'r941_min_high_g360',
    'r941_min_sup_g507',
    'r941_prom_hac_g507',
    'r941_prom_high_g360',
    'r941_prom_sup_g507'
]


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Specific options for the hybrid assembly pipeline.
    """
    ploidy: int | None = dataclasses.field(default=1, metadata={'choices': [1, 2]})
    ont_qual: str | None = dataclasses.field(default='nano-corr', metadata={'choices': ['nano-corr', 'nano-hq', 'nano-raw']})
    flye_meta: bool | None = dataclasses.field(default=None, metadata={'help': 'enables the --meta option for the Flye assembly'})
    expected_genome_size: str | None = dataclasses.field(default='', metadata={'help': 'Expected genome size'})
    ont_basecalling_model: str | None = dataclasses.field(default='r941_prom_sup_g507', metadata={'choices': BASECALLING_MODELS})
    filtlong_keep_percent: int | None = dataclasses.field(default=None)
    unicycler: bool = dataclasses.field(default=False)
    clair3_long_indel: bool = dataclasses.field(default=None)
    sniffles_min_svlen: int | None = dataclasses.field(default=None)
    sniffles_mapq: int  | None= dataclasses.field(default=None)
    sniffles_min_support: int  | None= dataclasses.field(default=None)
    freebayes_min_alternate_fraction: float | None = dataclasses.field(default=0.5)
    freebayes_min_alternate_count: int | None = dataclasses.field(default=10)


class MainHybridAssemblyPipeline(BasePipe):
    """
    Main class to run the hybrid assembly pipeline.
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
            name="Hybrid assembly pipeline",
            version="1.0",
            script_in=in_,
            script_out=out,
            opts=opts,
            snakefile=SNAKEFILE_MAIN,
        )
        self._opts_custom = opts_custom
        self._data_models = pd.read_table(TSV_BASECALLING_MODELS)

    def run(self) -> None:
        """
        Runs the hybrid assembly pipeline.
        :return: None
        """
        path_config = self.__create_snakemake_config_data()
        self.run_snakefile(path_config)

    def __create_snakemake_config_data(self) -> str:
        """
        Creates a Snakemake configuration file.
        :input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Config file data
        """
        config_data = self.get_config_data()

        # Assembly steps
        config_data['assembly_steps'] = ['Flye', 'Medaka', 'Polypolish', 'Pypolca']
        config_data['base_assemblies'] = ['flye']
        if self._opts_custom.unicycler:
            config_data['assembly_steps'].extend([
                'Unicycler', 'Medaka-Unicycler', 'Polypolish-Unicycler', 'Pypolca-Unicycler'])
            config_data['base_assemblies'].extend(['unicycler'])

        # Template data
        with open(CONFIG_DATA) as handle_in:
            basepipeutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    clair3_haploid_precise=True if self._opts_custom.ploidy == 1 else False,
                    clair3_long_indel=True if self._opts_custom.clair3_long_indel is not None else False,
                    clair3_model_path=str(self._get_clair3_model(self._opts_custom.ont_basecalling_model)),
                    freebayes_min_alternate_count=self._opts_custom.freebayes_min_alternate_count,
                    freebayes_min_alternate_fraction=self._opts_custom.freebayes_min_alternate_fraction,
                    freebayes_ploidy=self._opts_custom.ploidy,
                    expected_size=self._opts_custom.expected_genome_size if self._opts_custom.expected_genome_size else False,
                    is_meta=True if self._opts_custom.flye_meta else False,
                    medaka_model=self._opts_custom.ont_basecalling_model,
                    nano_corr=True if self._opts_custom.ont_qual == 'nano-corr' else False,
                    nano_hq=True if self._opts_custom.ont_qual == 'nano-hq' else False,
                    nano_raw=True if self._opts_custom.ont_qual == 'nano-raw' else False,
                    sniffles_mapq=self._opts_custom.sniffles_mapq if self._opts_custom.sniffles_mapq is not None else 40,
                    sniffles_min_support=self._opts_custom.sniffles_min_support if self._opts_custom.sniffles_min_support is not None else 'auto',
                    sniffles_min_svlen=self._opts_custom.sniffles_min_svlen if self._opts_custom.sniffles_min_svlen is not None else 35
                )))
        return snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

    def _get_clair3_model(self, medaka_model: str) -> Path:
        """
        Gets the most suitable clair3 model for the given medaka model.
        From documentation: the medaka model with the highest version equal to or less than the guppy version should
        be selected.
        :param medaka_model: Medaka basecalling model
        :return: Path to the Clair3 model
        """
        clair3_model = next(
            r['clair3_model'] for r in self._data_models.to_dict('records') if r['medaka_model'] == medaka_model)
        logger.info(f"Best matching Clair3 model for medaka model '{medaka_model}': {clair3_model}")
        path_model = Path(config.dir_db, 'clair3', 'models', clair3_model)
        if not path_model.exists():
            raise ValueError(f'Clair3 model not found: {path_model}')
        return path_model


@click.command(name='hybrid_assembly_pipeline', short_help='Hybrid assembly pipeline (ONT + Illumina)')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Hybrid assembly pipeline.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(**cliutils.from_kwargs(Options, kwargs))
    pipeline = MainHybridAssemblyPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
