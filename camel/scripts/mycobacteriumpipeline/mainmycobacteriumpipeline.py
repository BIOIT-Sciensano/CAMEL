#!/usr/bin/env python
import dataclasses
import tempfile
from importlib.resources import files
from pathlib import Path

import click
import yaml

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.core.snakemake import snakemakeutils
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import logger, initialize_logging
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basepipe.basepipe import BasePipe
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.tools.pipelines.mycobacterium.bamaddcustomtag import BAMAddCustomTag
from camel.scripts.mycobacteriumpipeline import SNAKEFILE_MAIN
from camel.snakefiles import variant_calling

CUSTOM_ANALYSES = [
    '51snp',
    'amr',
    'cgmlst',
    'confindr',
    'csb_rd',
    'hsp65',
    'human_read_scrubbing'
    'kraken2',
    'mlst',
    'ncbi_16s',
    'rmlst',
    'snp_lineage',
    'snpit',
    'spoligotyping',
]

@dataclasses.dataclass(frozen=True)
class Options(model.BaseOptions):
    """
    Pipeline-specific options.
    """
    analyses: list[str] = dataclasses.field(default_factory=list)
    output_bam: Path | None = dataclasses.field(default=None, metadata={
        'help': 'Output path for the mapping to the reference genome (BAM)'})


class MainMycobacteriumPipeline(BasePipe):
    """
    Main class to run the Mycobacterium pipeline.
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
            name='Mycobacterium pipeline',
            title='<i>Mycobacterium</i> pipeline',
            version='1.3',
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
        with open(str(files('camel').joinpath('scripts/mycobacteriumpipeline/config_data.yml'))) as handle:
            yaml_text = handle.read()
        yaml_text = yaml_text.format(
            COV_MAX=self._script_opts.cov_max,
            QC_SCHEME='cgmlst' if 'cgmlst' in  self._opts_custom.analyses else 'mlst',
            EXPORT_BAM=self._script_opts.include_bam
        )
        data_template = yaml.safe_load(yaml_text)
        self._script_out.dir.mkdir(parents=True, exist_ok=True)

        # Add the base config data
        config_data = self.get_config_data()
        basepipeutils.dict_merge(config_data, data_template)
        config_data['analyses'] = self._opts_custom.analyses
        config_data['sequence_typing']['options'] = {'method': self._script_opts.typing_method}
        config_data['gene_detection']['options'] = {'method': self._script_opts.gene_detection_method}

        # Map to reference genome instead of assembly
        if self._script_in.type_ is model.InputType.ILLUMINA:
            config_data['quality_checks']['forced'] = ['map_rate_ref_illumina', 'cov_ref_illumina']
            config_data['quality_checks']['skipped'] = ['map_rate_assembly_illumina', 'cov_assembly_illumina']
        if self._script_in.type_ is model.InputType.ONT:
            config_data['quality_checks']['forced'] = ['map_rate_ref_ont', 'cov_ref_ont']
            config_data['quality_checks']['skipped'] = ['map_rate_assembly_ont', 'cov_assembly_ont']

        path_config = snakepipelineutils.generate_config_file(config_data, self._script_opts.working_dir)

        # Run the Snakefile
        self.run_snakefile(path_config)
        self._export_assembly()
        self._export_bam()

    def _export_bam(self) -> None:
        """
        Exports the BAM output to the specified output location (optional).
        :return: None
        """
        if self._opts_custom.output_bam is None:
            logger.debug('Not exporting BAM output')
            return
        path_io = variant_calling.get_bam({
            'input': {'type': self._script_in.type_.value},
            'working_dir': self._script_opts.working_dir
        })
        bam_input = snakemakeutils.load_object(Path(self._script_opts.working_dir, path_io))

        # Add a custom tag with the sample name for PACU
        with tempfile.TemporaryDirectory(prefix='camel_', dir=config.dir_temp) as dir_:
            add_tag = BAMAddCustomTag()
            add_tag.add_input_files({'BAM': bam_input})
            add_tag.update_parameters(
                output=str(self._opts_custom.output_bam), name='PACU_name', value=self._script_in.name)
            add_tag.run(Path(str(dir_)))
        logger.info(f'Output BAM file copied to: {self._opts_custom.output_bam}')

@click.command(name='mycobacterium_pipeline', short_help='Pipeline for the complete characterization of Mycobacterium tuberculosis complex isolates')
@basescriptutils.add_input_opts()
@basescriptutils.add_output_opts
@basescriptutils.add_general_opts
@click.option('--analyses', type=str, help=f"Comma-separated list of analyses to run ({', '.join(CUSTOM_ANALYSES)})")
@cliutils.add_click_options_from_dataclass(Options, skip=['analyses'])
def main(**kwargs) -> None:
    """
    Pipeline for the complete characterization of Mycobacterium tuberculosis complex (MTBC) isolates.
    """
    script_input = basescriptutils.parse_script_input(kwargs)
    script_out = basescriptutils.parse_script_output(kwargs)
    script_opts = basescriptutils.parse_script_opts(kwargs)
    custom_opts = Options(
        analyses=kwargs['analyses'].split(',') if kwargs['analyses'] else [],
        output_bam=kwargs['output_bam']
    )
    pipeline = MainMycobacteriumPipeline(script_input, script_out, script_opts, custom_opts)
    pipeline.run()


if __name__ == '__main__':
    initialize_logging()
    main()
