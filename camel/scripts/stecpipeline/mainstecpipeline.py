#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import yaml

from camel.app.scriptutils import mainscriptutils
from camel.app.scriptutils.reportpipeline import ReportPipeline
from camel.app.core.snakemake import snakepipelineutils
from camel.app.loggers import initialize_logging
from camel.scripts.stecpipeline import CONFIG_DATA, SNAKEFILE_MAIN


class MainSTECPipeline(ReportPipeline):
    """
    Main class to run the STEC pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'rmlst', 'amrfinder', 'resfinder4', 'ncbi_stress', 'mlst_pasteur', 'mlst_warwick',
        'cgmlst', 'plasmidfinder', 'mob_suite', 'serotype', 'virulencefinder', 'innuendo_cgmlst',
        'human_read_scrubbing', 'variant_calling']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('STEC pipeline', '1.2', SNAKEFILE_MAIN, args)

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_input_files()
        path_config = self.__construct_config_file(input_files)
        self._run_snakemake_main(str(path_config))
        self._export_assembly()

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> Path:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = [key for key in MainSTECPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    coverage_max=self._args.cov_max,
                    qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst_warwick',
                    export_bam='true' if self._args.report_include_bam else 'false',
            )))

        # Read trimming
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library
        return Path(snakepipelineutils.generate_config_file(config_data, self._args.working_dir))

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainSTECPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    initialize_logging()
    main = MainSTECPipeline()
    main.run()
