import abc
from typing import Any

from camel.app.core.command import Command


class BaseDependencyService(metaclass=abc.ABCMeta):
    """
    Baseclass for dependency service.
    """

    def __init__(self) -> None:
        """
        Initializes the service.
        """
        pass

    @abc.abstractmethod
    def setup_environment(self, tool_data: dict[str, Any]) -> None:
        """
        Setup an environment.
        :param tool_data: Tool data
        :return: None
        """
        pass

    @abc.abstractmethod
    def load_environment(self, command: Command, tool_data: dict[str, Any]) -> str:
        """
        Loads an environment.
        :param command: Command to run
        :param tool_data: Tool data
        :return: Command with environment loaded
        """
        pass

    @abc.abstractmethod
    def is_available(self, tool_data: dict[str, Any]) -> bool:
        """
        Checks if the target environment is available.
        :param tool_data: Tool data
        :return: True if available, False otherwise
        """
        pass
