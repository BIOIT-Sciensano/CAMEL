import os
import re
import shutil
import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest
import yaml

from camel.app.config import config
from camel.app.core.reports import reportutils
from camel.app.loggers import logger


@pytest.fixture
def running_dir() -> Generator[Path]:
    """
    Creates a temporary working directory for tests.
    Automatically removed unless CAMEL_KEEP_TEST_DIRS=1.
    """
    dir_test = Path(tempfile.mkdtemp(prefix='camel_', dir=config.dir_temp))
    logger.debug(f"Directory for testing: {dir_test}")
    yield dir_test

    if os.environ.get("CAMEL_KEEP_TEST_DIRS") == '1':
        logger.debug("Keeping working directory (CAMEL_KEEP_TEST_DIRS)")
        return

    if dir_test.exists():
        logger.debug(f"Removing working directory: {dir_test}")
        shutil.rmtree(dir_test, ignore_errors=True)


def get_test_file_dir(*args) -> Path:
    """
    Retrieves the directory with test files.
    :return: Path to test file directory
    """
    dir_test = Path(config.dir_testdata, *args)
    if not dir_test.exists() or not dir_test.is_dir():
        raise FileNotFoundError(f"Cannot find test file directory: {dir_test}")
    return dir_test


def get_reference_file_dir(*args) -> Path:
    """
    Retrieves the directory with reference files.
    """
    dir_ref = Path(config.dir_db, "refgenomes", *args)
    if not dir_ref.exists() or not dir_ref.is_dir():
        raise FileNotFoundError(f"Cannot find reference file directory: {dir_ref}")
    return dir_ref


def verify_output_files(tool, key: str, nb_files: int = 1) -> None:
    """
    Verifies if the specified number of output files with the given key are created correctly.
    :param tool: Tool to check
    :param key: Output key
    :param nb_files: Number of expected output files
    :return: None
    """
    assert key in tool.tool_outputs, f"Key '{key}' missing from tool outputs"

    actual = len(tool.tool_outputs[key])
    assert actual == nb_files, f"Unexpected number of tool outputs found: {actual}, expected {nb_files}"

    for i, output in enumerate(tool.tool_outputs[key]):
        output_file_path = output.path
        assert output_file_path.exists(), f"Output file '{key}' (index: {i}) does not exist"
        assert output_file_path.stat().st_size > 0, f"Output file '{key}' (index: {i}) is empty"


def export_report_section(report_section, dir_out: Path) -> Path:
    """
    Exports the report section to the output directory.
    """
    path_out = dir_out / "report_section.html"
    path_out.parent.mkdir(exist_ok=True, parents=True)
    report = reportutils.init_report(
        path_out=path_out,
        key="report_section",
        title="Exported report section",
        dir_out=dir_out)
    report.add_html_object(report_section)
    report.save()
    logger.info(f"Report section exported to: {path_out}")
    return path_out

def extract_from_yaml(path_in: Path, key: str, placeholders: dict[str, str] | None = None) -> dict:
    """
    Extracts the given top-level section from the input YAML file.
    :param path_in: Input YAML file
    :param key: Section key
    :param placeholders: Placeholders to replace in the YAML file
    :return: Parsed section data
    """
    with open(path_in) as handle:
        lines = handle.read().splitlines()

    try:
        start = next(i for i, l in enumerate(lines) if l.startswith(f'{key}:'))
    except StopIteration:
        raise ValueError(f"Key '{key}' not found in YAML file: {path_in}")

    # End is the next top-level key after start
    pattern = re.compile(r'^[A-Za-z0-9_]+:')
    for i in range(start + 1, len(lines)):
        if pattern.match(lines[i]):
            end = i
            break
    else:
        end = len(lines)

    yaml_region = '\n'.join(lines[start:end])
    yaml_region = yaml_region.format(**placeholders if placeholders is not None else {})
    return yaml.safe_load(yaml_region)[key]
