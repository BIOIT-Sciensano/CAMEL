import abc
from typing import TypeVar, Generic

from camel.app.loggers import logger
from camel.app.scriptutils import model

TInput = TypeVar("TInput", bound=model.BaseInput)
TOutput = TypeVar("TOutput", bound=model.BaseOutput)
TOpts = TypeVar("TOpts", bound=model.BaseOptions)


class BaseScript(Generic[TInput, TOutput, TOpts], metaclass=abc.ABCMeta):
    """
    Base class for main scripts.
    """

    def __init__(
        self,
        name: str,
        version: str,
        script_in: TInput,
        script_out: TOutput,
        script_opts: TOpts,
        title: str | None = None,
    ) -> None:
        """
        Initializes the main script.
        :param name: Script name
        :param version: Script version
        :return: None
        """
        self._name: str = name
        self._version: str = version
        self._title: str = title
        self._script_in: TInput = script_in
        self._script_out: TOutput = script_out
        self._script_opts: TOpts = script_opts

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def title(self) -> str:
        """
        Returns the script title (can contain additional formatting).
        :return: Title
        """
        return self._title if self._title else self.name

    @property
    def version(self) -> str:
        """
        Returns the pipeline version.
        :return: Version
        """
        return self._version

    def info(self) -> dict[str, str]:
        """
        Returns the script information as a dictionary.
        :return: Dictionary with script information
        """
        return {
            'name': self.name,
            'title': self.title,
            'version': self.version,
        }

    def log_run(self) -> None:
        """
        Logs the tool execution.
        :return: None
        """
        logger.info(f"Running script: {self.name} (v{self.version})")

    @abc.abstractmethod
    def _execute(self) -> None:
        """
        Executes the script.
        This method should be overwritten by the subclasses.
        :return: None
        """
        raise NotImplementedError("Should be implemented by subclasses")

    def run(self) -> None:
        """
        Runs the script.
        """
        self.log_run()
        self._script_in.validate()
        self._execute()
