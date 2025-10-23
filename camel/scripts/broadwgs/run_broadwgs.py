#!/usr/bin/env python
import argparse
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import yaml

from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.pipeline.pipeline import Pipeline
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.core.snakemake import snakepipelineutils
from camel.scripts.broadwgs import (
    CONFIG_DATA,
    REFERENCES,
    SLURM_SUBMIT,
    SNAKEFILE_MAIN,
    TOOL_DATA,
)


class MainBroadWGSPipeline:
    """
    Main class to run the GATK best practices workflow
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        self._name = 'Broad institute WGS GATK Best Practices'
        self._version = '1.0'
        self._snakefile = SNAKEFILE_MAIN
        self._args = MainBroadWGSPipeline._parse_arguments(args)
        self._working_dir = Path(self._args.working_dir)
        self._final_output_dir = self._working_dir / "output"
        self._pipeline = Pipeline(self._name, Camel.get_instance(), self._args.log, self._args.log)

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def version(self) -> str:
        """
        Returns the pipeline version.
        :return: Version
        """
        return self._version

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        self._process_cmdline_args()

        if self._args.config is not None:
            config_file = self._args.config
        else:
            config_file = self.__construct_config_file()
        with open(config_file) as handle_in:
            config_data = yaml.load(handle_in.read(), Loader=yaml.SafeLoader)

        try:
            if self._args.slurm:
                logger.info("Running pipeline with Slurm: resources and threads input parameters ignored.")
                threads = 1
                resources = None
                slurm_args = {'cluster': self._args.slurm}
                if self._args.slurm == f'python3 {SLURM_SUBMIT} {{dependencies}}':
                   slurm_args.update(config_data["slurm"]["slurm_args"])
                logger.debug(slurm_args)
                SnakePipelineUtils.run_snakemake(self._snakefile, config_file, [], self._working_dir, threads,
                                                 resources, slurm_args)
            else:
                resources = dict([arg.split(",") for arg in self._args.resources])
                SnakePipelineUtils.run_snakemake(self._snakefile, config_file, [], self._working_dir, self._args.threads, resources)

            log_file = self._working_dir / 'camel.log'
            if log_file.exists():
                shutil.copyfile(str(log_file), str(Path(self._final_output_dir) / 'camel.log'))
            logger.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            if self._pipeline.keep_error_log:
                self._pipeline.log_error_to_file(err)
            raise err

    def _process_cmdline_args(self) -> None:
        """
        Input files: prepare ToolIOFiles, symlink
        :return: None
        """
        base_path = self._args.input

        # Create directory
        dir_links = self._working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True, exist_ok=True)

        # Determine link names + generate ToolIOFiles
        links = []
        for file in base_path:
            #add extension
            r1 = Path(f"{file}R1.fastq.gz")
            r2 = Path(f"{file}R2.fastq.gz")
            #append to links for symlinking
            links.append(r1)
            links.append(r2)
            #ToolIOFile
            io_files = [ToolIOFile(r1), ToolIOFile(r2)]
            SnakemakeUtils.dump_object(io_files, dir_links / f"{Path(file).name}.fastq.gz.io")

        # Link files
        for path_orig in links:
            link_name = path_orig.name
            path_new = dir_links / link_name
            logger.debug(f"Symlinking input file: {link_name}")
            if path_new.is_symlink():
                path_new.unlink()
            path_new.symlink_to(path_orig)

    def __construct_config_file(self) -> str:
        """
        Constructs the configuration file
        :return: Configuration file
        """
        # pipeline information and cmd line arguments

        config_data = {
            'pipeline': {
                'name': self._name,
                'version': f"{self._version}",
            },
            'sample': self._args.sample,
            'input_basenames': [Path(f).name for f in self._args.input],
            'working_dir': str(self._working_dir),
            'final_output_dir': str(self._final_output_dir),
            'debug': self._args.debug,
            'no_qc': self._args.no_qc
        }

        # add data from config yml
        if self._args.config_data is not None:
            config_data_file = Path(self._args.config_data)
        else:
            config_data_file = CONFIG_DATA
        with config_data_file.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # add snakemake parameters from command line
        config_data["params_smk"]['resources'] = self._args.resources
        config_data["params_smk"]['threads'] = self._args.threads

        # add data from user specified references yml if defined
        if self._args.references is not None:
            reference_file = Path(self._args.references)
        else:
            reference_file = REFERENCES
        with reference_file.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # add tool parameters from tool_data yml
        if self._args.tool_data is not None:
            tool_data = Path(self._args.tool_data)
        else:
            tool_data = TOOL_DATA
        with tool_data.open() as handle_in:
            config_data.update(yaml.load(handle_in.read(), Loader=yaml.SafeLoader))

        # Update read length
        config_data["rule_params"]["qc"]["picard_wgs_metrics"]["read_length"] = self._args.read_length
        config_data["rule_params"]["qc"]["picard_raw_wgs_metrics"]["read_length"] = self._args.read_length

        return snakepipelineutils.generate_config_file(config_data, self._working_dir)

    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser(description='Run the GATK best practices pipeline.')

        # input data
        parser.add_argument('-i', dest = "input", nargs="+", help="FastQ input files")
        parser.add_argument('--working-dir', dest = "working_dir", default = str(Path('.').absolute()), type=str, help='Working directory')
        parser.add_argument('--sample', type = str, help = 'Sample name')
        parser.add_argument('--read-length', dest = "read_length", type = int, default = 250, help = 'Read length. Default: 250')

        # input files
        parser.add_argument('--config', type = str, help = "Path to full config yml. Used to rerun analysis using an existing config file")
        parser.add_argument('--config-data', dest = "config_data", type=str, help="Path to config data yml")
        parser.add_argument('--references', type = str, help = "Path to references yml")
        parser.add_argument('--intervals', type = str, help = 'Path to interval files directory')
        parser.add_argument('--tool-data', dest = "tool_data", type = str, help = "Path to tool data yml")

        # logs
        parser.add_argument('--log', action='store_true', help="If this flag is set, config file and error logs are kept")

        # Snakemake parameters
        parser.add_argument('--threads', type = int, default = 4, help = "Snakemake parameter: number of cores")
        parser.add_argument('--resources', type = str, default = "mem_mb,5000", nargs = "+", help =
        "Snakemake parameter: resources. Key-value pair separated by a comma, e.g. mem_mb,1000. Multiple pairs allowed")

        # other
        parser.add_argument('--debug', dest = "debug", action = 'store_true')
        parser.add_argument('--slurm', dest = "slurm", nargs = "?", const = f'python3 {SLURM_SUBMIT} {{dependencies}}', type = str)
        parser.add_argument('--no-qc', dest = "no_qc", action = 'store_true')
        parser.set_defaults(debug = False)

        return parser.parse_args(args)


if __name__ == '__main__':
    initialize_logging()
    main = MainBroadWGSPipeline()
    main.run()
