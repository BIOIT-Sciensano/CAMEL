from importlib.resources import files
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel


class CamelConfig(BaseModel):
    """
    Config for the CAMEL project.
    """
    dir_testdata: Optional[Path] = None
    dir_temp: Optional[Path] = None
    dir_db: Path
    dir_error_logs: Path | None
    timestamp_format: str = '%Y-%d-%m %H:%M:%S'
    dependency_service: str = 'lmod'
    logging_fmt: str = '%(asctime)s - %(module)15s - %(levelname)7s - %(message)s'
    date_fmt: str = '%d/%m/%Y - %X'
    # Logging
    dir_configs: Optional[Path] = None
    dir_logs: Optional[Path] = None


def load_yaml_config(config_path: Path) -> dict:
    """
    Loads the configuration from a YAML file.
    :param config_path: Configuration path.
    :return: Config data as a dictionary.
    """
    with config_path.open('r') as file:
        return yaml.safe_load(file)


config_path_ = files('camel').joinpath('config/camel.yml')
config_yaml = yaml.safe_load(config_path_.read_text())
config = CamelConfig(**config_yaml)
