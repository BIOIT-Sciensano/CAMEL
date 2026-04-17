#!/usr/bin/env python
import dataclasses
from importlib.resources import files

import click
import yaml

from camel.app.cli import cliutils
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
from camel.scripts.bacilluspipeline import SNAKEFILE_MAIN

CUSTOM_ANALYSES = {
    "common": [
        "rmlst",
        "plasmidfinder",
        "mob_suite",
        "vfdb_core",
        "amrfinder",
        "kraken2",
        "confindr",
        "human_read_scrubbing",
        "variant_calling",
    ],
    "cereus": ["btyper", "mlst_cereus", "cgmlst_cereus"],
    "subtilis": ["fastani", "mlst_subtilis", "gmo", "straingst"],
}


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """

    species: str = dataclasses.field(metadata={'choices': list(CUSTOM_ANALYSES.keys())})
    analyses: list[str] = dataclasses.field(default_factory=list)
    mobsuite_contig_report: bool = dataclasses.field(
        default=False, metadata={'help': 'Export contig report for MOB-suite'}
    )


class MainBacillusPipeline(BasePipe):
    """
    Main class to run the Bacillus pipeline.
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
            name='Bacillus pipeline',
            title='<i>Bacillus</i> pipeline',
            version='1.0',
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
        with open(
            str(files('camel').joinpath('scripts/bacilluspipeline/config_data.yml'))
        ) as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            DB_ROOT=config.dir_db,
            EXPORT_BAM=self._script_opts.include_bam,
            MOBSUITE_CONTIG_REPORT=self._opts_custom.mobsuite_contig_report,
            SPECIES=self._opts_custom.species,
            QC_SCHEME=self.__get_qc_typing_scheme(),
        )
        data_template = yaml.safe_load(yaml_text)

        # Add the base config data
        config_data = self.get_config_data()
        config_data['species_selected'] = self._opts_custom.species

        # Add the analyses
        config_data["analyses_selected"] = []
        for key in self._opts_custom.analyses:
            if (
                key not in CUSTOM_ANALYSES["common"]
                and key not in CUSTOM_ANALYSES[self._opts_custom.species]
            ):
                logger.warning(
                    f"Analysis '{key}' not supported for species '{self._opts_custom.species}'"
                )
                continue
            config_data["analyses_selected"].append(key)

        # Update the config data
        basepipeutils.dict_merge(config_data, data_template)
        config_data['sequence_typing']['options'] = {
            'method': self._script_opts.typing_method
        }
        config_data['gene_detection']['options'] = {
            'method': self._script_opts.gene_detection_method
        }

        # Resolve species-specific values
        config_data = basepipeutils.resolve_config(config_data, self._opts_custom.species)
        return config_data

    def _validate_config_data(self, config_data: dict) -> bool:
        """
        Validates the config data.
        :param config_data: Config data
        :return: True if valid, False otherwise
        """
        self.check_dbs(config_data)
        return True

    def __get_qc_typing_scheme(self) -> str:
        """
        Returns the typing scheme used for QC.
        :return: Typing scheme
        """
        if self._opts_custom.species == "cereus":
            return (
                "cgmlst_cereus"
                if 'cgmlst_cereus' in self._opts_custom.analyses
                else "mlst_cereus"
            )
        return "mlst_subtilis"

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
    name='bacillus_pipeline',
    short_help='Pipeline for the complete characterization of Bacillus isolates',
)
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option(
    '--analyses',
    type=str,
    help=f"Comma-separated list of analyses to run ({', '.join(a for _, analyses in CUSTOM_ANALYSES.items() for a in analyses)})",
)
@cliutils.add_click_options_from_dataclass(Options, skip=['analyses'])
def main(**kwargs) -> None:
    """
    Pipeline for the complete characterization of Bacillus isolates.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
        species=kwargs['species'],
    )
    pipeline = MainBacillusPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.prepare_input()
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
