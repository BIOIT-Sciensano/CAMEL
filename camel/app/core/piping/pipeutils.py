import tempfile
from dataclasses import dataclass
from pathlib import Path

from camelcore.app.command import Command

from camel.app.config import config
from camel.app.core.piping.toolpipeable import ToolPipeable


@dataclass
class PipedTool:
    """
    Holder class for piped tools.
    """
    tool: ToolPipeable
    command: str
    stderr_path: Path


def _combine_dependencies(tools: list[ToolPipeable]) -> str:
    """
    Combines the dependencies for tools.
    :param tools: List of tools.
    :return: Combined lmod command
    """
    # Combine dependencies
    versions_by_name = {}
    for tool in tools:
        for dependency in tool.dependencies:
            name, version = dependency.split('/')
            if name not in versions_by_name:
                versions_by_name[name] = []
            if version not in versions_by_name[name]:
                versions_by_name[name].append(version)

    # Check if there are conflicts
    for name, versions in sorted(versions_by_name.items()):
        if len(versions) > 1:
            raise RuntimeError(f"Incompatible dependencies: {name} [{', '.join(versions)}]")

    # Create novel module load command
    return f"module load {' '.join([f'{name}/{versions[0]}' for name, versions in versions_by_name.items()])}"


def run_as_pipe(tools: list[ToolPipeable], dir_: Path) -> list[PipedTool]:
    """
    Runs a set of tools as a single pipe.
    :param tools: List of tools
    :param dir_: Working directory
    :return: List of piped tools
    """
    # Collect separate commands
    piped_tools = []
    for i, tool in enumerate(tools):
        stderr_path = Path(tempfile.NamedTemporaryFile(dir=config.dir_temp, prefix='stderr_', suffix='.txt').name)
        piped_tools.append(PipedTool(
            tool,
            f'{tool.prepare_pipe(dir_, i > 0, i != len(tools) - 1).command} 2> {stderr_path}',
            stderr_path
        ))

    # Construct and run full command
    lmod_command = _combine_dependencies(tools)
    piped_command = ' | '.join([piped_tool.command for piped_tool in piped_tools])
    full_command = Command(f'set -eo pipefail; {lmod_command}; {piped_command}')
    full_command.run(dir_)
    if not full_command.exit_code == 0:
        raise ValueError("error executing pipe")

    # Run after pipe steps
    for i, piped_tool in enumerate(piped_tools):
        with piped_tool.stderr_path.open() as handle:
            stderr = handle.read()
        piped_tool.tool.process_pipe(stderr, i == len(tools) - 1)
        piped_tool.tool.informs['_command'] = full_command.command

    return piped_tools
