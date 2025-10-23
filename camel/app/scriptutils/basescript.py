import abc
import argparse
from pathlib import Path
from collections.abc import Sequence


class BaseScript(metaclass=abc.ABCMeta):
    """
    Base class for main scripts.
    """

    def __init__(self, name: str, version: str, snakefile: Path | None) -> None:
        """
        Initializes the main script.
        :param name: Script name
        :param version: Script version
        :param snakefile: Snakefile path
        :return: None
        """
        self._name = name
        self._version = version
        self._snakefile = snakefile
        self._args: argparse.Namespace | None = None

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments(args: Sequence[str] | None) -> argparse.Namespace:
        """
        Parses the command line arguments. Should be implemented by the subclasses.
        :param args: Arguments (optional)
        :return: None
        """
        pass

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
