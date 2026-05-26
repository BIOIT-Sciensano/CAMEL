import json
import pickle
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.io.tooliovalue import ToolIOValue

from camel.app.core.tool import Tool
from camel.app.loggers import logger


class IOEncoder(json.JSONEncoder):
    """
    This class is used to encode IO objects in JSON format.
    """

    def default(self, obj: Any) -> dict:
        """
        Default encoding function.
        :param obj: Object to be encoded
        :return: Encoded object as a dictionary
        """
        if isinstance(obj, ToolIOFile):
            return {
                'io_type': 'ToolIOFile',
                'path': str(obj.path),
                'logged': obj.is_logged,
            }
        elif isinstance(obj, ToolIODirectory):
            return {
                'io_type': 'ToolIODirectory',
                'path': str(obj.path),
                'logged': obj.is_logged,
            }
        elif isinstance(obj, ToolIOValue):
            return {
                'io_type': 'ToolIOValue',
                'value': obj.value,
                'logged': obj.is_logged,
            }
        return super().default(obj)


def io_hook(obj: dict) -> object | dict:
    """
    Method to decode IO objects from JSON.
    :param obj: Dictionary with object properties
    :return: Decoded object
    """
    if 'io_type' not in obj:
        return obj
    if obj['io_type'] == 'ToolIOFile':
        return ToolIOFile(path=Path(obj['path']), logged=obj['logged'])
    elif obj['io_type'] == 'ToolIODirectory':
        return ToolIODirectory(path=Path(obj['path']), logged=obj['logged'])
    elif obj['io_type'] == 'ToolIOValue':
        return ToolIOValue(value=obj['value'], logged=obj['logged'])
    else:
        raise ValueError(f"Unsupported IO type: {obj['io_type']}")


def dump_object(obj: Any, path: Path) -> None:
    """
    Dumps an object in a pickle.
    :param obj: Object to dump
    :param path: Path to store the pickle
    :return: None
    """
    logger.debug(f"Dumping object '{obj}' in file '{path}'")
    if path.suffix == '.io':
        with open(path, 'w') as handle:
            json.dump(obj, handle, cls=IOEncoder, indent=2)
    elif path.suffix == '.iob':
        with path.open('wb') as handle:
            # noinspection PyTypeChecker
            pickle.dump(obj, handle)
    else:
        raise ValueError('Unsupported extension, expected .io (JSON) or .iob (pickle)')


def load_object(path: Path) -> Any:
    """
    Loads the object from the given IO file.
    :param path: Path to IO file
    :return: Object
    """
    logger.debug(f"Loading object from file '{path}'")
    if path.name.endswith('.io'):
        with path.open() as handle:
            return json.load(handle, object_hook=io_hook)
    elif path.name.endswith('.iob'):
        with path.open('rb') as handle:
            return pickle.load(handle)
    raise ValueError(f'Cannot load IO object from: {path}: invalid extension')


def add_io_input(tool: Tool, key: str, path: Path, optional: bool = False) -> None:
    """
    Adds an IO input (pickle or json) to a tool. For optional input whose value is empty, it is skipped.
    :param tool: Tool
    :param key: Key
    :param path: IO path (.io for json, .iob for pickle)
    :param optional: True for optional input, False otherwise
    :return: None
    """
    logger.debug(
        f"Adding IO input with key '{key}' from file '{path}' to tool '{tool.name}'"
    )
    value = load_object(path)
    if optional and len(value) == 0:
        logger.debug(f"Optional Input '{key}' empty, skipped")
    else:
        tool.add_input_files({key: value})


def add_io_inputs(
    tool: Tool,
    snake_input: Any,
    keys: list[str] | None = None,
    excluded_keys: list[str] | None = None,
    optionals: list[str] | None = None,
) -> None:
    """
    Adds pickled inputs from the snakemake input. If 'optionals' is specified, any optional input in that
    list will be skipped if its value is empty (no input file).
    :param tool: Tool
    :param snake_input: Snakemake input
    :param keys: Keys to add. If None, all keys are added
    :param excluded_keys: For the keys in this list the files are not added
    :param optionals: list of keys specifying optional inputs
    :return: None
    """
    logger.info("Adding pickled inputs from snakemake input")

    # Set variables
    optionals = [] if optionals is None else optionals
    keys = snake_input.keys() if keys is None else keys

    # Add keys
    for key in keys:
        logger.debug(f"Adding input '{key}'")
        if key not in snake_input.keys():
            raise KeyError(f"Key '{key}' not found in snakemake input")

        # Check if key was excluded
        if (excluded_keys is not None) and (key in excluded_keys):
            continue

        # Load the value
        if key in optionals and len(snake_input[key]) == 0:
            logger.debug(f"Optional Input '{key}' empty, skipped")
            continue

        value = load_object(Path(snake_input[key]))
        if key.startswith('INFORMS'):
            inform_key = '_'.join(key.split('_')[1:])
            tool.add_input_informs({inform_key: value})
            logger.debug(f"Informs '{value}' added")
        else:
            tool.add_input_files({key: value})
            logger.debug(f"Input '{value}' added")


def dump_io_output(tool: Tool, key: str, path: Path) -> None:
    """
    Dumps a tool output to an IO file.
    :param tool: Tool
    :param key: Key
    :param path: IO path
    :return: None
    """
    logger.debug(
        f"Dumping output with key '{key}' from tool '{tool}' to Camel IO pickle '{path}'"
    )
    if key not in tool.tool_outputs:
        raise KeyError(f"Tool '{tool.name}' has no output '{key}'")
    dump_object(tool.tool_outputs[key], path)


def dump_io_outputs(
    tool: Tool,
    snake_output: Any,
    keys: list[str] | None = None,
    ignore_missing_output: bool = False,
) -> None:
    """
    Dumps the tool outputs in pickles.
    :param tool: Tool
    :param snake_output: Snake output
    :param keys: Keys to dump
    :param ignore_missing_output: If False, an error is raised when an output is not generated
    :return: None
    """
    logger.info("Dumping tool outputs")
    if keys is None:
        keys = snake_output.keys()
    for key in keys:
        if key in tool.tool_outputs:
            dump_object(tool.tool_outputs[key], Path(snake_output[key]))
        elif key == 'INFORMS':
            dump_object(tool.informs, Path(snake_output[key]))
        else:
            message = f"Output '{key}' not generated"
            if ignore_missing_output is True:
                logger.warning(message)
            else:
                raise ValueError(message)


def update_param_if_not_none(
    tool: Tool, key: str, params: dict | Any, transform: Callable = None
) -> None:
    """
    Updates tool parameters if the value is not None.
    :param tool: Tool instance
    :param key: Parameter key
    :param params: Parameter dictionary
    :param transform: Transformation function
    :return: None
    """
    value = params.get(key)
    if value is None:
        return
    tool.update_parameters(**{key: value if transform is None else transform(value)})


def export_to_tsv(data: list[tuple[str, str | int | float]], output: Path):
    """
    Exports the data to a TSV file.
    :param data: list of tuples, where each tuple represents a key-value pair to be included in the TSV
    :param output: path to the output TSV file
    """
    with output.open('w') as handle:
        for key, value in data:
            handle.write(f'{key}\t{value}')
            handle.write('\n')


def export_to_json(
    data: list[tuple[str, str | int | float]],
    output: Path,
    main_key: str | None = None,
) -> None:
    """
    Exports the data to a JSON file.
    :param data: A list of tuples, each representing a key-value pair to be included in the JSON output.
    :param output: Path to the output JSON file.
    :param main_key: The primary key under which the data is organized in the JSON file; if none is provided, the
        data is written directly into the file.
    :return: None
    """
    # Convert numpy integers to regular ints
    for idx, (k, v) in enumerate(data):
        if not isinstance(v, np.int64):
            continue
        data[idx] = (k, int(v))
    json_content = {main_key: dict(data)} if main_key else dict(data)
    with output.open('w') as handle:
        json.dump(json_content, handle, indent=2)


def export_summary(
    data: list[tuple[str, str | int | float]],
    path_out: Path,
    ext: str,
    json_main_key: str | None = None,
) -> None:
    """
    Exports the summary output in the given format.
    :param data: Data to export
    :param path_out: Output path
    :param ext: Format (.tsv / .json)
    :param json_main_key: Main key to use for the JSON output
    :return: None
    """
    if ext == 'json':
        export_to_json(data, path_out, json_main_key)
    elif ext == 'tsv':
        export_to_tsv(data, path_out)
    else:
        raise ValueError(f'Invalid ext: {ext}')


def convert_list_to_dict(
    data: list[list[str]], headers: list[str]
) -> list[dict[str, str]]:
    """
    Converts the input list of lists to a list of dicts and adds headers.
    :param data: list of lists
    :param headers: headers that need to be added
    :return: list of dicts where each dict contains the headers
    """
    df = pd.DataFrame(data, columns=headers)
    df.replace(np.nan, '-', inplace=True)
    # noinspection PyTypeChecker
    return df.to_dict(orient='records')


def sanitize_numpy(obj: Any) -> Any:
    """
    Recursively convert numpy scalars in nested structures to Python scalars.
    This method can be used for sanitizing informs before storing then in JSON format.
    :return: None
    """
    # Convert numpy scalars (np.int64, np.float64, np.bool_, etc.)
    if isinstance(obj, np.generic):
        return obj.item()

    # Dict → sanitize keys & values
    elif isinstance(obj, dict):
        return {sanitize_numpy(k): sanitize_numpy(v) for k, v in obj.items()}

    # List → sanitize each element
    elif isinstance(obj, list):
        return [sanitize_numpy(v) for v in obj]

    # Tuple → sanitize then rebuild tuple
    elif isinstance(obj, tuple):
        return tuple(sanitize_numpy(v) for v in obj)

    # Set → sanitize then rebuild set
    elif isinstance(obj, set):
        return {sanitize_numpy(v) for v in obj}

    # Everything else stays the same
    return obj


def get_rule_dir(output) -> Path:
    """
    Derives the working directory for a rule from its output paths.
    All outputs must share the same parent directory.
    :param output: Snakemake output Namedlist
    :return: Parent directory of the output files
    :raises ValueError: If output is empty or outputs span more than one directory
    """
    dirs = {Path(p).parent for p in output}
    if not dirs:
        raise ValueError('Cannot derive working directory from empty output.')
    if len(dirs) > 1:
        raise ValueError(
            f'Rule outputs span multiple directories: {dirs}. '
            f'All outputs must share a single working directory.'
        )
    return dirs.pop()
