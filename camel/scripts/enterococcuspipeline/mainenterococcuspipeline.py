#!/usr/bin/env python
import dataclasses
from pathlib import Path
from typing import Any

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

CUSTOM_ANALYSES = [
    'amrfinder',
    'bacmet',
    'cgmlst',
    'confindr',
    'human_read_scrubbing',
    'kraken2',
    'lrefinder',
    'mlst',
    'mob_suite',
    'plasmidfinder',
    'resfinder4',
    'rmlst',
    'variant_calling'
    'vfdb_core',
    'virulencefinder',
]

DATA_BY_SPECIES = {
    'faecalis': {
        'amrfinder_species': 'Enterococcus_faecalis',
        'cgmlst_db': str(Path(config.dir_db, 'sequence_typing/enterococcus_faecalis/cgmlst')),
        'full_name': 'Enterococcus faecalis',
        'reference': {
            'fasta': str(Path(config.dir_db, 'refgenomes/Enterococcus_faecalis/KB944666.1.fasta')),
            'gff3': str(Path(config.dir_db, 'refgenomes/Enterococcus_faecalis/KB944666.1.gff3')),
            'name': 'KB944666.1',
            'url': 'https://www.ncbi.nlm.nih.gov/nuccore/KB944666.1',
        },
        'mlst_db': str(Path(config.dir_db, 'sequence_typing/enterococcus_faecalis/mlst')),
        'resfinder4_species': 'Enterococcus faecalis',
    },
    'faecium': {
        'amrfinder_species': 'Enterococcus_faecium',
        'cgmlst_db': str(Path(config.dir_db, 'sequence_typing/enterococcus_faecium/cgmlst')),
        'full_name': 'Enterococcus faecium',
        'reference': {
            'fasta': str(Path(config.dir_db, 'refgenomes/Enterococcus_faecium/CP038996.1.fasta')),
            'gff3': str(Path(config.dir_db, 'refgenomes/Enterococcus_faecium/CP038996.1.gff3')),
            'name': 'CP038996.1',
            'url': 'https://www.ncbi.nlm.nih.gov/nuccore/CP038996.1/',
        },
        'mlst_db': str(Path(config.dir_db, 'sequence_typing/enterococcus_faecium/mlst')),
        'resfinder4_species': 'Enterococcus faecium',
    },
    'spp': {
        'amrfinder_species': None,
        'disabled_assays': ['mlst', 'cgmlst', 'variant_calling'],
        'full_name': 'Enterococcus spp.',
        'reference': {
            'gc': 37.4,
            'size': 2_973_380
        },
        'resfinder4_species': None,
    },
}


@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """
    species: str = dataclasses.field(metadata={'choices': list(DATA_BY_SPECIES.keys())})
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
            name='Enterococcus pipeline',
            title='<i>Enterococcus</i> pipeline',
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
            AMRFINDER_SPECIES=DATA_BY_SPECIES[self._opts_custom.species]['amrfinder_species'],
            CGMLST_DB=DATA_BY_SPECIES[self._opts_custom.species].get('cgmlst_db'),
            COV_MAX=self._script_opts.cov_max,
            EXPORT_BAM='true' if self._script_opts.include_bam else 'false',
            K2_LEVEL='S' if self._opts_custom.species != 'spp' else 'G',
            K2_NAME=DATA_BY_SPECIES[self._opts_custom.species][
                'full_name'] if self._opts_custom.species != 'spp' else 'Enterococcus',
            MLST_DB=DATA_BY_SPECIES[self._opts_custom.species].get('mlst_db'),
            QC_SCHEME='cgmlst' if 'cgmlst' in self._opts_custom.analyses else 'rmlst',
            REF_FASTA=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('fasta', 'null'),
            REF_GC=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('gc', 'null'),
            REF_GFF3=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('gff3', 'null'),
            REF_NAME=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('name', 'null'),
            REF_SIZE=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('size', 'null'),
            REF_URL=DATA_BY_SPECIES[self._opts_custom.species].get('reference', {}).get('url', 'null'),
            RESFINDER4_SPECIES=DATA_BY_SPECIES[self._opts_custom.species]['resfinder4_species'],
        )
        data_template = yaml.safe_load(yaml_text)
        self._script_out.dir.mkdir(parents=True, exist_ok=True)

        # Add the base config data
        config_data = self.get_config_data()
        basepipeutils.dict_merge(config_data, data_template)
        config_data['analyses'] = self._opts_custom.analyses
        config_data['species'] = self._opts_custom.species
        config_data['sequence_typing']['options'] = {'method': self._script_opts.detection_method}
        config_data['gene_detection']['options'] = {'method': self._script_opts.detection_method}

        # Additional MLST scheme for E. faecium
        if (self._opts_custom.species == 'faecium') and ('mlst' in self._opts_custom.analyses):
            config_data['analyses'].append('mlst_bezdicek')

        # Disable species-specific assays for generic Enterococcus
        config_data['is_generic'] = self._opts_custom.species == 'spp'
        if self._opts_custom.species == 'spp':
            self._update_config_for_generic_spp(config_data)
        config_data['selected_species'] = DATA_BY_SPECIES[self._opts_custom.species]['full_name']

        # Create the config file
        path_config = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run the Snakefile
        snakepipelineutils.run_snakemake(
            snakefile=self._snakefile,
            config_path=path_config,
            targets=[],
            working_dir=self._script_opts.working_dir,
            threads=self._script_opts.threads)
        self._export_assembly()

    def _update_config_for_generic_spp(self, config_data: dict[str, Any]) -> None:
        """
        Updates the config file with specific adaptation for generic enterococcus.
        :param config_data: Configuration data
        :return: None
        """
        # Disable incompatible assays
        disabled_assays = DATA_BY_SPECIES['spp']['disabled_assays']
        config_data['analyses'] = [a for a in config_data['analyses'] if a not in disabled_assays]
        logger.warning(f"Generic 'Enterococcus' selected as species, disabling assays: {', '.join(disabled_assays)}")

        # Disable species specific AMR detection
        config_data['amrfinder']['species'] = None
        config_data['resfinder4']['species'] = None
        config_data['resfinder4']['point'] = False

        # Change the typing scheme for the QC check (no cgMLST is available)
        logger.warning("cgMLST is not available for generic 'Enterococcus', using rMLST for the QC check.")
        config_data['quality_checks']['typing_scheme'] = 'rmlst'


@click.command(name='enterococcus_pipeline', short_help='Pipeline for the complete characterization of Enterococcus isolates')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option('--analyses', type=str, help=f"Comma-separated list of analyses to run ({', '.join(CUSTOM_ANALYSES)})")
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
    pipe_script = MainEnterococcusPipeline(script_input, script_out, script_opts, custom_opts)
    pipe_script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
