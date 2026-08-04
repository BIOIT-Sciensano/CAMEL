import logging
from pathlib import Path

from camel.app.config import config

logger = logging.getLogger('camel')

_VERBOSE_LOGGERS = ('camel', 'camelcore')


def initialize_logging() -> None:
    """
    Initializes the logging.
    :return: None
    """
    root = logging.getLogger()

    # Avoid attaching the handlers more than once
    if any(getattr(h, 'name', None) == 'console' for h in root.handlers):
        return
    formatter = logging.Formatter(config.logging_fmt)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.DEBUG)
    console_handler.name = 'console'
    root.addHandler(console_handler)

    # File handler (camel.log file)
    file_handler = logging.FileHandler('camel.log')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    file_handler.name = 'file'
    root.addHandler(file_handler)

    # Only set the camel-related loggers to DEBUG level
    root.setLevel(logging.WARNING)
    for name in _VERBOSE_LOGGERS:
        logging.getLogger(name).setLevel(logging.DEBUG)


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
    # Attached to the root logger so the step log captures both CAMEL and camel-core records.
    logging.getLogger().addHandler(handler)


def detach_step_handlers() -> None:
    """
    Detaches & closes the handlers used for this step.
    """
    root = logging.getLogger()
    to_remove = []
    for handler in root.handlers:
        if (handler.name is None) or (not handler.name.startswith('step_')):
            continue
        to_remove.append(handler)
        handler.close()
    for handler in to_remove:
        root.removeHandler(handler)


initialize_logging()
