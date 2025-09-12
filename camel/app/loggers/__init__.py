import logging
from pathlib import Path

from camel.app.config import config

logger = logging.getLogger('camel')


def initialize_logging() -> None:
    """
    Initializes the logging.
    :return: None
    """
    if logger.hasHandlers():
        return
    formatter = logging.Formatter(config.logging_fmt)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.DEBUG)
    console_handler.name = 'console'
    logger.addHandler(console_handler)

    # File handler (camel.log file)
    file_handler = logging.FileHandler('camel.log')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    # General logging level
    logger.setLevel(logging.DEBUG)
    logger.propagate = False


def attach_step_handler(dir_out: Path, level: int) -> None:
    """
    Attaches a FileHandler for the given step.
    :param dir_out: Output directory (for storing logs)
    :param level: Logging level
    :return: None
    """
    dir_out.mkdir(exist_ok=True, parents=False)
    path_out = dir_out / f'{logging.getLevelName(level).lower()}.log'
    handler = logging.FileHandler(path_out)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(config.logging_fmt))
    handler.name = f'step_{level}'
    logger.addHandler(handler)


def detach_step_handlers() -> None:
    """
    Detaches & closes the handlers used for this step.
    """
    to_remove = []
    for handler in logger.handlers:
        if (handler.name is None) or (not handler.name.startswith('step_')):
            continue
        to_remove.append(handler)
        handler.close()
    for handler in to_remove:
        logger.removeHandler(handler)


initialize_logging()
