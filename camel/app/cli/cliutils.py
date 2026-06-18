import dataclasses
import importlib
import pkgutil
import traceback
import typing
from collections.abc import Callable
from pathlib import Path
from types import UnionType
from typing import Any, get_args, get_origin

import click
from click.testing import Result

from camel.app.loggers import logger


def list_scripts(script_dir: Path) -> typing.Iterator[tuple[Path, str]]:
    """
    Lists the available scripts in the given directory.
    :param script_dir: Directory containing the scripts
    :return: Iterator over directory and script names
    """
    for dir_ in sorted(script_dir.iterdir()):
        if not dir_.is_dir() or dir_.name.startswith("_"):
            continue
        for _, name, _ in pkgutil.iter_modules([str(dir_)]):
            if not name.startswith("main"):
                continue
            yield dir_, name


def load_script_module(script_name: str, script_dir: Path) -> Any:
    """
    Loads the given script module.
    :param script_name: Script name
    :param script_dir: Script directory
    :return: Imported module
    """
    module_path = f"camel.scripts.{script_dir.name}.{script_name}"
    return importlib.import_module(module_path)


def type_to_click(field_type: Any) -> Any:
    """
    Maps a dataclass type to a click type.
    :param field_type: Dataclass type
    :return: Click type
    """
    origin = get_origin(field_type)
    if origin is list:
        inner_type = get_args(field_type)[0]
        if inner_type == Path:
            return click.Path(path_type=Path)
        return inner_type
    if field_type == Path:
        return click.Path(path_type=Path)
    if field_type == bool:
        # boolean -> is_flag options should be used
        return None
    return field_type


def unwrap_optional(field_type) -> tuple[Any, bool]:
    """
    Unwraps an optional type by removing the 'None' type.
    :param field_type: Optional type
    :return: Unwrapped type
    """
    origin = get_origin(field_type)
    optional = False
    if origin is UnionType:
        optional = type(None) in get_args(field_type)
        args = [t for t in get_args(field_type) if t is not type(None)]
        if args:
            return args[0], optional
    return field_type, optional


def add_click_options_from_dataclass(
    dataclass_type: Any, skip: list[str] | None = None
) -> Callable:
    """
    Adds click options for all fields in the given dataclass.
    :param dataclass_type: Input dataclass
    :param skip: Fields to skip
    :return: Decorator to add click options.
    """

    def decorator(f: Callable) -> Callable:
        """
        Returns the decorator that adds options for all fields.
        :param f: Input function
        :return: Decorator
        """
        for field in reversed(dataclasses.fields(dataclass_type)):
            # Check if the field should be skipped
            if skip is not None and field.name in skip:
                continue

            # Unwrap optional types and map to click types
            field_type, optional = unwrap_optional(field.type)
            click_type = type_to_click(field_type)

            # List with possible values
            if field.metadata.get("choices") is not None:
                click_type = click.Choice(field.metadata["choices"])

            # Explicitly defined option type
            if field.metadata.get("type") is not None:
                click_type = field.metadata["type"]

            # Help text
            help_text = field.metadata.get("help", "")
            option_args = {
                "default": (
                    field.default if field.default != dataclasses.MISSING else None
                ),
                "help": help_text,
                "required": not optional and (field.default is dataclasses.MISSING),
                "show_default": field.metadata.get('show_default', False),
            }

            # Boolean flags
            if field.type == bool:
                option_args["is_flag"] = True
                option_args.pop("default", None)
                option_args.pop("required", None)

            # List types
            elif get_origin(field.type) is list:
                option_args["multiple"] = True

            # Add the option
            f = click.option(
                f"--{field.name.replace('_','-')}", type=click_type, **option_args
            )(f)
        return f

    return decorator


def from_kwargs(cls: Any, kwargs: dict, skip: list[str] | None = None) -> dict:
    """
    Extracts the fields from the given kwargs.
    :param cls: Dataclass type
    :param kwargs: Keyword arguments
    :param skip: List of fields to skip
    :return: Filtered kwargs
    """
    field_names = {f.name for f in dataclasses.fields(cls)}
    if skip is not None:
        field_names -= set(skip)
    for k, v in kwargs.items():
        if isinstance(v, Path):
            kwargs[k] = v.expanduser().resolve().absolute()
    return {k: v for k, v in kwargs.items() if k in field_names}


def invoke(script: Callable, args: list[str]) -> Result:
    """
    Invokes the target script with Click.
    :param script: Script to run
    :param args: Script arguments
    :return: Result
    """
    runner = click.testing.CliRunner()
    # noinspection PyTypeChecker
    result = runner.invoke(script, args, catch_exceptions=False)
    if result.exception is not None:
        logger.warning("".join(traceback.format_exception(result.exception)))
    logger.info(f"stdout: {result.stdout}")
    logger.info(f"stderr: {result.stdout}")
    return result
