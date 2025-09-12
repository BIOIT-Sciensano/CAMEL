#!/usr/bin/env python
import argparse
from importlib.resources import files
from pathlib import Path

from camel.app.loggers import initialize_logging
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from collections.abc import Sequence

class MainDummy:
    """
    Class to run the dummy script.
    """

    def __init__(self, args: Sequence[str] | None = None) -> None:
        """
        Initializes the class.
        :param args: Arguments (optional)
        :return: None
        """
        self._args = MainDummy._parse_args(args)

    @staticmethod
    def _parse_args(args: Sequence[str] | None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments to parse (optional)
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument('--threads', type=int, default=1)
        return parser.parse_args(args)

    def run(self) -> None:
        """
        Runs this script.
        :return: None
        """
        path_config = SnakePipelineUtils.generate_config_file({
            'input': str(self._args.fasta),
            'output': str(self._args.output)
        }, Path('.'))
        path_snakefile = Path(str(files('camel').joinpath('scripts/dummy/dummy.smk')))
        SnakePipelineUtils.run_snakemake(
            path_snakefile, path_config, [self._args.output], Path('.'), self._args.threads
        )

def run(args: Sequence[str] | None = None) -> None:
    """
    Runs the pipeline.
    :param args: Command line arguments
    :return: None
    """
    script = MainDummy(args)
    script.run()


if __name__ == '__main__':
    initialize_logging()
    run()
