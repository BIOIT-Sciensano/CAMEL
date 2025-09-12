from pathlib import Path
from typing import Union

import yaml

from camel.app.loggers import initialize_logging, logger
from camel.config import MAIN_CONFIG


class Camel:
    """
    Main class for camel.
    """

    _current_instance = None
    _logger_is_initialized = False

    def __init__(self, tool_parameter_loc: str = None) -> None:
        """
        Initializes a CAMEL system.
        :param tool_parameter_loc: Location of tool parameter YAML files
        """
        initialize_logging()
        with open(MAIN_CONFIG) as f:
            self._config = yaml.safe_load(f)

        if self._config.get('tool_service', 'db') == 'yaml' and 'tool_parameter_loc' not in self._config:
            self._config['tool_parameter_loc'] = tool_parameter_loc

        commit_hash = Camel.get_commit_hash()
        logger.debug(f"CAMEL commit hash: {commit_hash if commit_hash is not None else 'Not available'}")

    @property
    def config(self) -> dict:
        """
        Returns the main config as specified in app/config/main.yml
        :return: Dict
        """
        return self._config

    @staticmethod
    def get_instance() -> 'Camel':
        """
        This method can be used to avoid recreating the CAMEL object multiple times.
        :return: Initialized CAMEL instance
        """
        if Camel._current_instance is None:
            Camel._current_instance = Camel()
        return Camel._current_instance

    @staticmethod
    def get_commit_hash() -> Union[str, None]:
        """
        Checks the commit hash of the CAMEL repository (if available).
        """
        path_version_txt = Path(__file__).parents[2] / 'VERSION'
        if not path_version_txt.exists():
            return None
        with path_version_txt.open() as handle:
            return handle.readline().strip()
