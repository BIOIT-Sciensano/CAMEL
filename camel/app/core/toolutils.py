import hashlib
import inspect
import json
import os
import shutil
import stat
from pathlib import Path
from typing import Any

from camelcore.app.command import Command
from camelcore.app.io.tooliodirectory import ToolIODirectory
from camelcore.app.utils import fileutils

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.loggers import logger


def check_input(tool: 'Tool', keys_required: list[str], keys_allowed: list[str] | None = None) -> None:  # noqa: F821
    """
    Checks if the provided tool input is valid.
    :param tool: Tool instance
    :param keys_required: Input keys that are required
    :param keys_allowed: Input keys that are allowed
    :return: None
    """
    # Check required keys
    for key in keys_required:
        if key not in tool.tool_inputs:
            raise InvalidToolInputError(f"'{key}' input is required")

    # Check allowed keys
    for key in tool.tool_inputs:
        if (keys_allowed is not None) and (key not in keys_allowed):
            raise InvalidToolInputError(f"'{key}' input is not allowed")


def check_tool_execution(tool: 'Tool', command: Command, exit_code: int = 0):  # noqa: F821
    """
    Checks if the tool executed successfully.
    :param tool: Tool instance
    :param command: Tool command
    :param exit_code: Expected exit code
    """
    if command.exit_code != exit_code:
        raise ToolExecutionError(tool.name, f"Error executing '{tool.name}', exit code: {command.exit_code}")


def show_parameters(tool) -> str:
    """
    Returns an overview of the current parameters as a string.
    :return: Parameter overview
    """
    parts = []
    for param_key, param in tool.params.items():
        if param.flag:
            parts.append(f'{param_key}: True')
        else:
            parts.append(f'{param_key}: {param.value}')
    return ', '.join(parts)


def get_test_cache_dir() -> Path | None:
    """
    Returns the test cache root directory if CAMEL_TEST_CACHE_DIR is set, None otherwise.
    :return: Cache root directory, or None if caching is disabled
    """
    dir_cache = os.environ.get('CAMEL_TEST_CACHE_DIR')
    return Path(dir_cache) if dir_cache else None


def snapshot_dir(dir_: Path) -> dict[Path, float]:
    """
    Snapshots the (relative path -> mtime) of all files in a directory. Taken before a command runs so that
    save_to_cache can capture only the files the command created or changed, which keeps caching correct when tools
    share a working directory.
    :param dir_: Directory to snapshot
    :return: Mapping of relative path to modification time
    """
    if not dir_.is_dir():
        return {}
    return {f.relative_to(dir_): f.stat().st_mtime for f in dir_.rglob('*') if f.is_file()}


def _normalize_param_value(value: Any, folder: Path) -> str:
    """
    Returns a stable string representation of a parameter value for cache-key hashing. Paths inside the tool's working
    directory are made relative, preventing per-run directories from making identical invocations produce different
    cache keys.
    :param value: Parameter value
    :param folder: Tool working directory
    :return: Normalized string representation
    """
    text = str(value)
    folder_text = str(folder)
    if folder_text and text.startswith(folder_text):
        return '<WORKDIR>' + text[len(folder_text) :]
    return text


def _directory_fingerprint(path: Path) -> str:
    """
    Returns a cheap directory fingerprint based on file metadata (name, size, mtime) rather than content. Symlinks are
    resolved first, so retargeting a link (e.g. latest) changes the fingerprint.
    :param path: Directory path
    :return: Fingerprint string
    """
    real_path = path.resolve()
    entries = sorted(
        (str(f.relative_to(real_path)), f.stat().st_size, f.stat().st_mtime)
        for f in real_path.rglob('*')
        if f.is_file()
    )
    return f'{real_path}:{entries}'


def _input_fingerprint(item: 'ToolIO') -> str:  # noqa: F821
    """
    Returns a fingerprint for a tool input used to build the cache key. Directories use a cheap metadata-based
    fingerprint, while other inputs use their own memoized hash.
    :param item: Tool input
    :return: Fingerprint string
    """
    if isinstance(item, ToolIODirectory):
        return _directory_fingerprint(item.path)
    return item.hash


def compute_cache_key(tool: 'Tool') -> str:  # noqa: F821
    """
    Computes a cache key for the current tool invocation from the tool's identity, its implementation source, its active
    parameters and a fingerprint of its inputs.
    :param tool: Tool instance
    :return: Cache key (SHA-256 hex digest)
    """
    source_path = Path(inspect.getfile(tool.__class__))
    data = {
        'tool': tool.name_full,
        'source_hash': fileutils.hash_file(source_path),
        'yaml_hash': fileutils.hash_file(tool.get_tool_data_path()),
        'params': {
            name: _normalize_param_value(param.value, tool.folder) for name, param in sorted(tool.params.items())
        },
        'inputs': {
            key: [_input_fingerprint(item) for item in items] for key, items in sorted(tool.tool_inputs.items())
        },
    }
    encoded = json.dumps(data, sort_keys=True).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def load_from_cache(command: Command, dir_: Path, cache_root: Path, tool_name: str, cache_key: str) -> bool:
    """
    Attempts to restore a previous command execution from the test cache by repopulating the working directory and
    restoring stdout, stderr, and the exit code. This allows downstream logic to behave as if the command had actually
    run. Tool informs are not restored, as cached values may reference the original working directory.
    :param command: Command that would have been executed
    :param dir_: Working directory
    :param cache_root: Test cache root directory
    :param tool_name: Tool name
    :param cache_key: Cache key for this invocation
    :return: True if the cache was hit, False otherwise
    """
    entry_dir = cache_root / fileutils.make_valid(tool_name) / cache_key[:2] / cache_key
    if not (entry_dir / '_SUCCESS').is_file():
        return False

    with (entry_dir / 'meta.json').open() as handle:
        meta = json.load(handle)

    workdir_cache = entry_dir / 'workdir'
    if workdir_cache.is_dir():
        for cached_file in workdir_cache.rglob('*'):
            if not cached_file.is_file():
                continue
            dest = dir_ / cached_file.relative_to(workdir_cache)
            shutil.copy2(cached_file, dest)
            # Update permissions so the file is writable by owner (cache is read only)
            dest.chmod(dest.stat().st_mode | stat.S_IWUSR)

    command._stdout = meta['stdout']
    command._stderr = meta['stderr']
    command._return_code = meta['exit_code']
    logger.info(f"Test cache hit for '{tool_name}' ({cache_key[:12]})")
    return True


def save_to_cache(
    command: Command, dir_: Path, cache_root: Path, tool_name: str, cache_key: str, files_before: dict[Path, float]
) -> None:
    """
    Stores the result of a real command execution into the test cache: the working directory files the command created
    or changed (relative to files_before) and the command's stdout, stderr and exit code. Only the delta is captured so
    that caching stays correct when tools share a working directory.
    :param command: Command that was executed
    :param dir_: Working directory
    :param cache_root: Test cache root directory
    :param tool_name: Tool name
    :param cache_key: Cache key for this invocation
    :param files_before: Snapshot of dir_ taken before the command ran (see snapshot_dir)
    :return: None
    """
    logger.debug(f"Saving to cache: '{tool_name}' ({cache_key})")
    entry_dir = cache_root / fileutils.make_valid(tool_name) / cache_key[:2] / cache_key
    if entry_dir.exists():
        shutil.rmtree(entry_dir)
    entry_dir.mkdir(parents=True)

    files_after = snapshot_dir(dir_)
    changed = [rel for rel, mtime in files_after.items() if files_before.get(rel) != mtime]
    for rel in changed:
        dest = entry_dir / 'workdir' / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        # Copy (not hardlink) so the cache entry is an independent inode from the live working file: a
        # hardlink shares content and permissions, so marking the entry read-only would also make the live
        # file read-only (breaking later in-place writes), and an in-place rewrite of the live file would
        # corrupt the cached copy. The independent copy is safe to mark read-only to protect the cache.
        shutil.copy2(dir_ / rel, dest)
        dest.chmod(0o444)

    meta = {'stdout': command.stdout, 'stderr': command.stderr, 'exit_code': command.exit_code}
    with (entry_dir / 'meta.json').open('w') as handle:
        json.dump(meta, handle)

    (entry_dir / '_SUCCESS').touch()
    logger.info(f"Stored test cache for '{tool_name}' ({cache_key[:12]}, {len(changed)} file(s))")
