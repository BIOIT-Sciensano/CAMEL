#!/usr/bin/env python
import importlib
import sys
import textwrap
from collections.abc import Sequence

from camel.app.loggers import logger, initialize_logging
from camel.version import __VERSION__


def main() -> None:
    """
    Main entry point to run CAMEL scripts or tools.
    :return: None
    """
    if len(sys.argv) < 2:  # noqa: PLR2004
        print(textwrap.dedent(
            """
            usage:
              camel
                run {SCRIPT}
                version
            """
        ))
        sys.exit(0)
    subcommand = sys.argv[1]
    match subcommand:
        case 'run':
            module = sys.argv[2]
            module_args: Sequence[str] = sys.argv[3:]
            try:
                mod = importlib.import_module(f"camel.scripts.{module}.main")
                mod.run(module_args)
            except ModuleNotFoundError:
                logger.warning(f"Error: Module '{module}' not found.")
                sys.exit(1)
        case 'version':
            print(f"CAMEL version: {__VERSION__}")


if __name__ == '__main__':
    initialize_logging()
    main()
