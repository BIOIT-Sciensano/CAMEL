import os
import subprocess
from pathlib import Path
from typing import Any

from camel.app.loggers import logger


class Command:
    """
    Helper class to run commands.
    """

    def __init__(self, command: str = None) -> None:
        """
        Initializes the command object.
        :param command: (optional) Command line call
        :return: None
        """
        self._stdout = None
        self._stderr = None
        self._procedure = None
        self._return_code = None
        self._command = command

    @property
    def stderr(self) -> str:
        """
        Returns the stderr from the command execution.
        :return: Standard error
        """
        return self._stderr

    @property
    def stdout(self) -> str:
        """
        Returns the stdout from the command execution.
        :return: Standard error
        """
        return self._stdout

    @property
    def returncode(self) -> int:
        """
        Returns the exit code from the command execution.
        :return: Exit code
        """
        return self._return_code

    @property
    def exit_code(self) -> int:
        """
        Returns the exit code from the command execution.
        """
        return self._return_code

    @property
    def command(self) -> str:
        """
        Returns the command line call.
        :return: Command line call
        """
        return self._command

    @command.setter
    def command(self, cmd: str) -> None:
        """
        Sets the command line call.
        :param cmd: Command
        :return: None
        """
        self._command = cmd

    def run(self, folder: Path, stderr_handle=subprocess.PIPE, disable_logging: bool = False, prefix: str | None = None,
            env: dict[str, Any] | None = None) -> None:
        """
        Runs the command given at command initialization
        :param folder: Folder where the command is executed
        :param stderr_handle: Handle for the standard error (e.g. PIPE or STDOUT)
        :param disable_logging: If True, logging is disabled
        :param prefix: If given, this prefix will be added to the commands
        :param env: Environment variables to be set for the command
        :return: None
        """
        command_str = self._command if prefix is None else f'{prefix}{self._command}'
        if disable_logging is False:
            logger.info(f'Executing command: {command_str}')
        if self.command is None:
            raise ValueError("Invalid command 'None'")
        self._procedure = subprocess.run(
            command_str,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            shell=True,
            executable='/bin/bash',
            env={**os.environ, **env} if env is not None else None,
            cwd=folder,
            text=True)
        self._stdout = self._procedure.stdout or ''
        self._stderr = self._procedure.stderr or ''
        self._return_code = self._procedure.returncode
        if disable_logging is False:
            logger.debug(f'stdout: {self._stdout}')
            logger.debug(f'stderr: {self._stderr}')

    def run_command(self, folder: Path, stderr_handle=subprocess.PIPE) -> None:
        """
        Runs the command given at command initialization
        :param folder: Folder where the command is executed
        :param stderr_handle: Handle for the standard error (e.g. PIPE or STDOUT)
        :return: None
        """
        import warnings
        warnings.warn(
            'The run_command method is deprecated, please use the run() method instead.', DeprecationWarning)
        self.run(folder, stderr_handle)
