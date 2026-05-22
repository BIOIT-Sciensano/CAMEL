from camelcore.app.command import Command

from camel.app.core.errors import InvalidToolInputError, ToolExecutionError


def check_input(tool: 'Tool', keys_required: list[str], keys_allowed: list[str] = None) -> None:  # noqa: F821
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


def check_tool_execution(tool: 'Tool', command: Command, exit_code: int=0):  # noqa: F821
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
