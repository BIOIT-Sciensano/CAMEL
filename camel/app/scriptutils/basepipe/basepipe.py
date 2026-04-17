import abc
import dataclasses
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, TypeVar

import click

from camel.app.core.utils import fastautils
from camel.app.dbs.dbutils import DBEntry
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseInput
from camel.app.core import errors
from camel.app.core.snakemake import snakepipelineutils, snakemakeutils
from camel.app.loggers import logger
from camel.snakefiles import assembly

TInput = TypeVar("TInput", bound=BaseInput)


class BasePipe(BaseScript[ScriptInput, ScriptOutput, ScriptOptions], metaclass=abc.ABCMeta):
    """
    Baseclass for pipelines that run a Snakefile.
    """

    def __init__(
        self,
        name: str,
        version: str,
        script_in: ScriptInput,
        script_out: ScriptOutput,
        opts: ScriptOptions,
        snakefile: Path | str,
        title: str | None = None,
    ):
        """
        Initializes this pipeline.
        :param name: Pipeline name
        :param version: Pipeline version
        :param script_in: Script input
        :param script_out: Script output
        :param snakefile: Snakefile to run
        :param title: Pipeline title
        :return: None
        """
        super().__init__(
            name=name,
            title=title,
            version=version,
            script_in=script_in,
            script_out=script_out,
            script_opts=opts,
        )
        self._snakefile = Path(snakefile)

    @staticmethod
    def check_analyses_option(value: str | None, allowed: list[str]) -> None:
        """
        Checks if the provided analyses option is valid.
        :param value: Analyses option string
        :param allowed: List of allowed analyses
        :return: None
        """
        if value is None:
            return
        for analysis in value.split(","):
            if analysis in allowed:
                continue
            raise click.UsageError(
                f"Invalid analysis option '{analysis}'. Allowed analyses: {', '.join(allowed)}"
            )

    def _validate_config_data(self, config_data: dict) -> bool:
        """
        Validates the config data.
        :param config_data: Config data
        :return: True if valid, False otherwise
        """
        self.check_dbs(config_data)
        return True

    def get_config_data(self) -> dict[str, Any]:
        """
        Returns the base config data.
        :return: Base config data
        """
        return {
            "input": self._script_in.to_dict(),
            "output": self._script_out.to_dict(),
            "script_info": self.info(),
            "working_dir": str(self._script_opts.working_dir),
            "read_trimming": {"method": self._script_opts.trimming_method},
        }

    def run_snakefile(self, path_config: Path | str) -> None:
        """
        Runs the pipeline snakefile.
        :param path_config: Path to the config file
        :return: None
        """
        # Clear existing Galaxy output files when HTML output is selected
        basepipeutils.prepare_galaxy_output(self._script_out.dir, self._script_out.html)

        # Path to the logfile
        log_file = self._script_opts.working_dir / "camel.log"
        try:
            snakepipelineutils.run_snakemake(
                self._snakefile,
                path_config,
                [],
                working_dir=self._script_opts.working_dir,
                threads=self._script_opts.threads,
            )
            logger.info("Pipeline finished successfully")
        except errors.SnakemakeExecutionError as err:
            if log_file.exists():
                path_log = basepipeutils.store_log_file(log_file, f'{self._name}_{self.version}', True)
                if path_log is not None:
                    message = f"Error executing Snakemake. Check log for more information: {path_log}"
                    raise RuntimeError(message)
            raise err

    def _export_assembly(self) -> None:
        """
        Exports the assembly to the specified output location (optional).
        :return: None
        """
        if self._script_out.fasta is None:
            logger.debug("Not exporting assembly")
            return
        path_io = self._script_opts.working_dir / assembly.OUTPUT_FASTA
        path_fasta = snakemakeutils.load_object(path_io)[0].path
        self._script_out.fasta.parent.mkdir(exist_ok=True, parents=True)
        fastautils.rename_sequences_regex(
            fasta_in=path_fasta,
            fasta_out=self._script_out.fasta,
            regex="",
            repl="",
            description=self._script_in.name,
        )
        logger.info(f"Output FASTA file exported to: {self._script_out.fasta}")

    def _execute(self) -> None:
        """
        Executes the pipeline.
        :return: None
        """
        basepipeutils.prepare_galaxy_output(self._script_out.dir, self._script_out.html)
        super()._execute()

    def check_dbs(self, data_template: dict) -> None:
        """
        Checks if the required databases are available.
        :param data_template: Template data
        :return: None
        """
        try:
            analyses = {key: data_template['analyses'][key] for key in data_template['analyses_selected']}
            logger.debug(f"Selected analyses: {', '.join(analyses.keys())}")
        except KeyError as err:
            raise ValueError(f"Analysis {err} is not defined for this pipeline")
        dbs = {key: DBEntry(**data) for key, data in data_template['dbs'].items()}

        # Check for unused DBs
        used_dbs = set(db for a in data_template['analyses'].values() for db in a["dbs"])
        defined_dbs = set(data_template["dbs"].keys())
        required_dbs = {key for key, db_info in data_template['dbs'].items() if db_info.get('required') is True}
        unused = defined_dbs - used_dbs - required_dbs
        if len(unused) > 0:
            logger.warning(f"The followings DBs are defined but not used: {', '.join(unused)}")

        # Check if the DBs for the specified analyses are available
        if not basescriptutils.check_dbs(dbs, analyses):
            logger.info("Required databases are missing, aborting pipeline")
            sys.exit(1)

    def prepare_input(self) -> None:
        """
        Prepares the script input by creating symlinks.
        :return: None
        """
        # Create the directory for the symlinks
        dir_links = self._script_opts.working_dir / 'input'
        dir_links.mkdir(parents=True, exist_ok=True)

        links_by_key = defaultdict(list)
        for key, path, link_name in self._script_in.get_symlinks():
            path_link_out = dir_links / link_name

            # Remove existing symlink
            if path_link_out.exists() or path_link_out.is_symlink():
                path_link_out.unlink()

            # Store the link
            path_link_out.symlink_to(path.resolve())
            links_by_key[key].append(path_link_out)

        # Convert lists back to single items (if needed)
        to_replace = {
            k: (v[0] if len(v) == 1 else v)
            for k, v in links_by_key.items()
        }
        self._script_in = dataclasses.replace(self._script_in, **to_replace)
