from pathlib import Path


def absolute_path_by_pathlib(path: str) -> Path:
    """
    Takes a relative path and returns the absolute path.
    :param path: Relative or absolute path
    :return: The resolved absolute Path object
    """
    return Path(path).expanduser().resolve()
