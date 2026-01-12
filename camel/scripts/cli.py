#! /usr/bin/env python3
import sys
from pathlib import Path

import click

from camel.app.cli import cliutils
from camel.app.config import config
from camel.app.dbs import dbutils
from camel.app.loggers import initialize_logging, logger
from camel.version import __VERSION__

SCRIPT_DIR = Path(__file__).parent

@click.group()
@click.version_option(__VERSION__, prog_name="CAMEL", message="%(prog)s %(version)s")
def cli() -> None:
    """
    Main command line interface.
    :return: None
    """
    pass

@cli.command(name='list')
def list_() -> None:
    """
    Lists the available scripts and pipelines.
    :return: None
    """
    for dir_, script_name in cliutils.list_scripts(SCRIPT_DIR):
        module = cliutils.load_script_module(script_name, dir_)
        cmd = getattr(module, "main", None)
        if isinstance(cmd, click.Command):
            click.echo(f"{cmd.name:<30} {cmd.short_help}")
        else:
            logger.warning(f"Script '{script_name}' does not define a Click command 'main'.")


@cli.group(context_settings=dict(max_content_width=120))
def run() -> None:
    """
    Runs a script or pipeline.
    :return: None
    """
    pass

def register_dynamic_subcommands() -> None:
    """
    Dynamically register subcommands for each script.
    :return: None
    """
    for dir_, script_name in cliutils.list_scripts(SCRIPT_DIR):
        module = cliutils.load_script_module(script_name, dir_)
        cmd = getattr(module, "main", None)
        if isinstance(cmd, click.Command):
            run.add_command(cmd)
        else:
            logger.warning(f"Script '{script_name}' does not define a Click command 'main'.")


register_dynamic_subcommands()

@cli.command(name='install')
@click.argument("yml", nargs=-1, type=click.Path(exists=True, path_type=Path))
def install(yml: list[Path]) -> None:
    """
    Installs tool dependencies.
    :param yml: Paths to YAML files to install.
    :return: None
    """
    logger.info("Installing is not supported yet")


@cli.command(name='download_dbs')
@click.pass_context
@click.option("--name", type=str, help='Pipeline name')
@click.option("--yml", type=click.Path(exists=True, path_type=Path), help='Pipeline YML file')
@click.option("--keys", type=str, help='List of keys to download (comma-separated), defaults to all')
@click.option("--force", is_flag=True, help='Force download, even if the DBs already exist')
@click.option("--out", type=click.Path(exists=True, path_type=Path), help='Output directory (defaults to value from config)')
@click.option("--threads", type=int, default=4, help='Number of threads to use for indexing')
def download_dbs(ctx, yml: Path, name: str, force: bool, keys: str | None, out: Path | None, threads: int) -> None:
    """
    Downloads databases for the CAMEL pipelines.
    """
    # Check if the options are valid
    if yml is None and name is None:
        raise click.UsageError("Either --yml or --name must be specified")

    # Determine the output directory
    dir_db = (out if out else config.dir_db)

    # Retrieve the YAML file path
    if name is not None:
        cmd_run = ctx.parent.command.commands.get('run')

        # Retrieve the pipeline command
        cmd_target = cmd_run.commands.get(name)
        if cmd_target is None:
            raise click.UsageError(f"Pipeline '{name}' not found")

        # Load the config data from the pipeline module
        module_name = cmd_target.callback.__module__
        module = sys.modules.get(module_name)
        try:
            path_yml = Path(getattr(module, 'CONFIG_DATA'))
        except AttributeError:
            raise click.UsageError(f"'{name}' does not define config data")
        if not path_yml.exists():
            raise FileNotFoundError(f"Config file '{path_yml}' not found")
    else:
        path_yml = yml

    # Install the DBs
    dbutils.download_dbs(
        yml=path_yml,
        force=force,
        keys=keys.split(',') if keys else None,
        dir_db=dir_db,
        threads=threads)


if __name__ == '__main__':
    initialize_logging()
    cli()
