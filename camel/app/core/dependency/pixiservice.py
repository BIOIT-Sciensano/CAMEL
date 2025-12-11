import hashlib
import shutil
from pathlib import Path
from typing import Any

from camel.app.config import config
from camel.app.core.command import Command
from camel.app.core.dependency.basedependencyservice import BaseDependencyService
from camel.app.core.errors import DependencyError
from camel.app.loggers import logger


class PixiService(BaseDependencyService):
    """
    Service for handling dependencies using Pixi.
    """

    def _get_dir_env(self, tool_data: dict[str, Any]) -> Path:
        """
        Returns the base directory for storing environments.
        :param tool_data: Tool data
        :return: Directory path
        """
        if config.dir_envs_pixi is None:
            raise ValueError("'dir_envs_pixi' is not set in the config")
        if 'conda' not in tool_data:
            raise ValueError("No 'conda' section found in tool data file")
        hash_str = hashlib.sha1(','.join(tool_data['conda']['packages']).encode()).hexdigest()[:8]
        return Path(config.dir_envs_pixi, f"env_{tool_data['conda']['name']}-{hash_str}")

    def setup_environment(self, tool_data: dict[str, Any]) -> None:
        """
        Setup an environment.
        :param tool_data: Tool data
        :return: None
        """
        # Create the directory
        dir_env = self._get_dir_env(tool_data)
        dir_env.mkdir(parents=True, exist_ok=True)

        # Create a new environment
        logger.info(f"Creating pixi environment: {tool_data['conda']['name']}")
        command = Command('pixi init --channel conda-forge --channel bioconda')
        command.run(dir_env)
        if not command.exit_code == 0:
            shutil.rmtree(dir_env)
            raise DependencyError(f'pixi init command failed: {command.stderr}')

        # Install packages
        pkgs = tool_data['conda']['packages']
        command = Command("pixi add {}".format(' '.join(f'"{p}"' for p in pkgs)))
        command.run(dir_env)
        if not command.exit_code == 0:
            shutil.rmtree(dir_env)
            raise DependencyError(f'pixi add command failed: {command.stderr}')

    def load_environment(self, command: Command, tool_data: dict[str, Any]) -> str:
        """
        Loads an environment.
        :param command: Command to run
        :param tool_data: Tool data
        :return: Command with environment loaded
        """
        logger.info('Loading environment using Pixi')
        dir_env = self._get_dir_env(tool_data)
        if '|' in command.command:
            command_sanitized = command.command.replace('"', r'\"')
            return f'pixi run --manifest-path {dir_env} bash -c "{command_sanitized}"'
        return f'pixi run --manifest-path {dir_env} {command.command}'

    def is_available(self, tool_data: dict[str, Any]) -> bool:
        """
        Checks if the target environment is available.
        :param tool_data: Tool data
        :return: True if available, False otherwise
        """
        try:
            dir_env = self._get_dir_env(tool_data)
        except ValueError:
            return False
        return dir_env.exists()
