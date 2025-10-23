from typing import Optional


class InvalidToolInputError(BaseException):
    """
    Error that is raised when the input to a tool is invalid.
    """
    pass

class ToolExecutionError(RuntimeError):
    """
    Error that is raised when the tool execution fails.
    """

    def __init__(self, tool_name: str, message: str) -> None:
        """
        Initializes the error.
        :param tool_name: Tool class in which the error occurred
        :param message: Message
        """
        super().__init__(f"Tool {tool_name} execution failed: {message}")

class InvalidParameterError(ValueError):
    """
    This is raised when an invalid parameter is supplied to a tool.
    """
    pass

class PipelineExecutionError(Exception):
    """
    Error that is raised when a pipeline cannot execute successfully.
    """
    pass


class SnakemakeExecutionError(RuntimeError):
    """
    This error is raised when a Snakemake execution error occurs.
    """

    def __init__(self, stdout: str, stderr: str, failed_rule: Optional[str] = None) -> None:
        """
        This class is raised when a snakemake error occurs.
        :param stdout: Standard output
        :param stderr: Error output
        """
        super().__init__(f"Failed at rule: {failed_rule if failed_rule else 'n/a'}")
        self.stdout = stdout
        self.stderr = stderr
        self.failed_rule = failed_rule
