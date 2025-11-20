#! /usr/bin/env python3
from pathlib import Path

import click

from camel.app.cli import cliutils
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


if __name__ == '__main__':
    initialize_logging()
    cli()
