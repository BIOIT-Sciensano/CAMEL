#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Any

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake import snakemakeutils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import human_read_scrubbing
from camel.scripts.ncbihumanreadscrubber import CONFIG_DATA
from camel.scripts.ncbihumanreadscrubber import SNAKEFILE_MAIN


class MainNcbiHumanReadScrubber(ReportPipeline):
    """
    Main class to run the NCBI human read scrubber tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('NCBI human read scrubber stand alone', '0.3', SNAKEFILE_MAIN, args)
        self._config_data: dict[str, Any] | None = None

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'NCBI human read scrubber'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._copy_output_files()

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = ['human_read_scrubbing']
        config_data['read_scrubbing'] = {}
        if self._args.export_removed_reads:
            config_data['read_scrubbing']['export_removed_reads'] = True
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        self._config_data = config_data
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    def _copy_output_files(self) -> None:
        """
        Copies the output files to the output directory.
        :return: None
        """
        # Copy the scrubbed reads
        output_files = snakemakeutils.load_object(self._args.working_dir / human_read_scrubbing.get_output_io(
            self._config_data))
        if self._args.input_type == 'illumina':
            shutil.copyfile(output_files[0].path, Path(self._args.output_dir / f'{self.sample_name}-scrubbed_R1.fastq.gz'))
            shutil.copyfile(output_files[1].path, Path(self._args.output_dir / f'{self.sample_name}-scrubbed_R2.fastq.gz'))
        elif self._args.input_type == 'ont':
            shutil.copyfile(output_files[0].path, Path(self._args.output_dir / f'{self.sample_name}-scrubbed.fastq.gz'))
        else:
            raise ValueError(f"Invalid input type: {self._args.input_type}")

        # Copy the removed reads
        if not self._args.export_removed_reads:
            return
        if self._args.input_type == 'ont':
            fq_removed = snakemakeutils.load_object(self._args.working_dir / human_read_scrubbing.get_removed('fastq_se'))
            if len(fq_removed) == 0:
                logger.warning('No removed reads found')
                return
            shutil.copyfile(fq_removed[0].path, self._args.output_dir / f'{self.sample_name}-removed.fastq.gz')
        elif self._args.input_type == 'illumina':
            fq_removed = snakemakeutils.load_object(self._args.working_dir / human_read_scrubbing.get_removed('fastq_pe'))
            if len(fq_removed) == 0:
                logger.warning('No removed reads found')
                return
            shutil.copyfile(fq_removed[0].path, self._args.output_dir / f'{self.sample_name}-removed_R1.fastq.gz')
            shutil.copyfile(fq_removed[1].path, self._args.output_dir / f'{self.sample_name}-removed_R2.fastq.gz')

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(argument_parser)
        argument_parser.add_argument(
            '--export-removed-reads', help="Export the removed reads", action='store_true')
        return argument_parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainNcbiHumanReadScrubber()
    main.run()
