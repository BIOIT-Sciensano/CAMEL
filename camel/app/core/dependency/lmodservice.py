from pathlib import Path
from typing import Any

from camel.app.core.command import Command
from camel.app.core.dependency.basedependencyservice import BaseDependencyService
from camel.app.loggers import logger


class LmodService(BaseDependencyService):
    """
    Service for loading tool dependencies using LMOD.
    """

    def setup_environment(self, tool_data: dict[str, Any]) -> None:
        """
        Setup an environment.
        :param tool_data: Tool data
        :return: None
        """
        raise NotImplementedError("Installing environments is not supported for LMOD.")

    def load_environment(self, command: Command, tool_data: dict[str, Any]) -> str:
        """
        Loads the LMOD environment.
        :param command: Command to run
        :param tool_data: Tool data
        :return: Command with environment loaded
        """
        logger.debug('Loading environment using LMOD')
        if len(tool_data.get('dependencies', [])) == 0:
            return command.command
        return ' '.join(['module load', *tool_data.get('dependencies', [])]) + '; '

    def is_available(self, tool_data: dict[str, Any]) -> bool:
        """
        Checks if the required dependencies are available.
        :param tool_data: Tool data for the target tool
        :return: None
        """
        for module in tool_data.get('dependencies', []):
            command = Command(f'module load {module}')
            command.run(Path.cwd(), disable_logging=True)
            if command.exit_code != 0:
                raise RuntimeError(f'Error loading module: {module}')
        return True
