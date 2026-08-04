import os
from importlib.resources import files
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel


class CamelConfig(BaseModel):
    """
    Config for the CAMEL project.
    """
    dir_testdata: Path | None = None
    dir_temp: Path | None = None
    dir_db: Path
    dir_error_logs: Path | None
    dir_envs_pixi: Path | None = None
    create_missing_envs: bool = False
    timestamp_format: str = '%Y-%d-%m %H:%M:%S'
    dependency_service: Literal['lmod', 'pixi'] = 'lmod'
    logging_fmt: str = '%(asctime)s - %(module)15s - %(levelname)7s - %(message)s'
    date_fmt: str = '%d/%m/%Y - %X'
    ftp_server: str | None = None
    # Logging
    dir_configs: Path | None = None
    dir_logs: Path | None = None


def load_yaml_config(config_path: Path) -> dict:
    """
    Loads the configuration from a YAML file.
    :param config_path: Configuration path.
    :return: Config data as a dictionary.
    """
    with config_path.open('r') as file:
        return yaml.safe_load(file)


def resolve_config_path() -> Path:
    """
    Resolves the path to the CAMEL configuration file.

    Looks in order at the ``$CAMEL_CONFIG`` environment variable, the user config
    directory (``$XDG_CONFIG_HOME`` or ``~/.config``)``/camel/main.yml``, and finally the
    ``main.yml`` shipped inside the package.
    :return: Path to the configuration file.
    """
    env_path = os.environ.get('CAMEL_CONFIG')
    if env_path:
        return Path(env_path)
    config_home = Path(os.environ.get('XDG_CONFIG_HOME') or Path.home() / '.config')
    user_path = config_home / 'camel' / 'main.yml'
    if user_path.is_file():
        return user_path
    packaged_path = Path(str(files('camel').joinpath('config/main.yml')))
    if packaged_path.is_file():
        return packaged_path
    raise FileNotFoundError(
        f"No CAMEL configuration found. Set $CAMEL_CONFIG or create {user_path} "
        "(copy camel/config/main.yml.sample as a starting point)."
    )


config_path_ = resolve_config_path()
config_yaml = load_yaml_config(config_path_)
config = CamelConfig(**config_yaml)
