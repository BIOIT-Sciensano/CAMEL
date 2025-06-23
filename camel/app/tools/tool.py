import abc
import inspect
from pathlib import Path
from typing import Optional, Union

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.error import InvalidToolInputError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.toolio import ToolIO
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.services.basetoolservice import BaseToolService
from camel.app.services.yamltoolservice import YAMLToolService


class Tool(metaclass=abc.ABCMeta):
    """
    Contains the common functionality of tools.
    """

    def __init__(self, name: str, version: str, camel: Optional[Camel] = None) -> None:
        """
        Initializes a tool.
        :param name: Tool name
        :param version: Tool version
        :param camel: CAMEL instance (optional)
        """
        logger.debug(f"Initializing tool: {name} {version}")
        self._name = name
        self._version = version
        self._tool_inputs: dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]] = {}
        self._tool_outputs: dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]] = {}
        self._informs = {'_name': self.name, '_version': self._version}
        self._input_informs = {}
        self._camel = camel
        self._tool_service = self.get_tool_service()
        self._tool_command = self._tool_service.get_tool_command()
        self._dependencies = self._tool_service.get_dependencies()
        self._parameters = self._tool_service.get_default_parameters()
        self._command = Command()
        self._folder = None

    @property
    def name(self) -> str:
        """
        Returns the name of this tool.
        :return: Name
        """
        return f'{self._name} {self._version}'

    @property
    def version(self) -> str:
        """
        Returns the tool version.
        :return: Version
        """
        return self._version

    @property
    def tool_id(self) -> int:
        """
        Returns the tool id.
        :return: Tool id
        """
        return self._tool_service.tool_id

    @property
    def tool_outputs(self) -> dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]]:
        """
        Returns the tool outputs.
        :return: Tool outputs
        """
        return self._tool_outputs

    @property
    def informs(self) -> dict:
        """
        Returns the tool informs.
        :return: Informs
        """
        return self._informs

    @property
    def dependencies(self) -> list[str]:
        """
        Returns a list of dependencies for this tool.
        :return: Dependencies
        """
        return self._dependencies

    @property
    def stdout(self) -> Optional[str]:
        """
        Returns the command line stdout.
        :return: Stdout
        """
        return self._command.stdout

    @property
    def stderr(self) -> Optional[str]:
        """
        Returns the command line stderr.
        :return: Stderr
        """
        return self._command.stderr

    @property
    def parameter_overview(self) -> str:
        """
        Returns an overview of the parameters as a string.
        :return: Parameters overview
        """
        return ', '.join([f"{p}: '{self._parameters[p].value}'" for p in sorted(self._parameters)]) if \
            len(self._parameters) > 0 else '/'

    @property
    def folder(self) -> Path:
        """
        Returns the folder the tool needs to run in.
        :return: Path to the running folder
        """
        return self._folder

    def add_input_files(self, input_files: dict[str, list[ToolIO]]) -> None:
        """
        Updates the input files for a tool.
        :param input_files: New input files
        :return: None
        """
        for key, items in input_files.items():
            if key in self._tool_inputs:
                self._tool_inputs[key] += items
            else:
                self._tool_inputs[key] = items

    def add_input_informs(self, informs: dict) -> None:
        """
        Updates the input informs for a tool.
        :param informs: New informs
        :return: None
        """
        self._input_informs.update(informs)

    def update_parameters(self, **kwargs: Union[str, int, None, bool, float, dict[str, Union[str, int, None, bool, float]]]) -> None:
        """
        Updates the parameters for this tool.
        :param kwargs: Parameters in key value format
        :return: None
        """
        for parameter_name, new_value in kwargs.items():
            parameter = self._tool_service.get_parameter(parameter_name)
            if not parameter:
                raise InvalidParameterError(f"{self._name} has no parameter '{parameter_name}'")
            if new_value is False:
                if parameter_name not in self._parameters:
                    logger.warning(f"Cannot disable parameter '{parameter_name}' (not present in parameters)")
                    continue
                logger.info(f"Disabling parameter: {parameter_name}")
                del(self._parameters[parameter_name])
            else:
                if new_value is True or new_value is None:
                    parameter.value = None
                else:
                    parameter.value = str(new_value)
                if parameter_name not in self._parameters:
                    logger.info(f"Parameter '{parameter_name}' added, value: {parameter.value}")
                else:
                    old_value = self._parameters[parameter_name].value
                    logger.info(f"Parameter '{parameter_name}' value '{old_value}' changed to '{new_value}'")
                self._parameters[parameter_name] = parameter

    def clear_parameters(self) -> None:
        """
        Clears all the parameters of the given tool.
        :return: None
        """
        logger.info(f"Removing {len(self._parameters)} parameters")
        self._parameters.clear()

    def run(self, folder: Path = Path.cwd()) -> None:
        """
        Runs this tool.
        :param folder: Folder to run the tool in.
        :return: None
        """
        self._folder = folder
        logger.info(f'Running tool {self.name}')
        logger.info(f'Working directory: {self._folder}')
        logger.info(f'Tool parameters: {self.parameter_overview}')
        self._check_parameters()
        self._check_input()
        self._execute_tool()
        self._check_output()

    def get_tool_data_path(self) -> Path:
        """
        Returns the path of the tool data for the tool with the given name and version.
        :return: Path
        """
        yaml_path = Path(inspect.getfile(self.__class__).replace('.py', '.yml'))
        if not yaml_path.is_file():
            raise FileNotFoundError(f"Tool data file for '{self.name}' not found ({yaml_path})")
        return yaml_path

    def get_tool_service(self) -> BaseToolService:
        """
        Returns the tool service for the tool with the given name and version.
        :return: Tool service
        """
        source = self._camel.config.get('tool_service', 'yaml')
        logger.debug(f'Retrieving tool service. Source = {source}')
        if source == 'db':
            raise DeprecationWarning("Parameter loading from database is deprecated, use YAML instead.")
        elif source == 'yaml':
            return YAMLToolService(self.get_tool_data_path())
        else:
            raise ValueError(f"Invalid 'tool_service' value: {source}")

    def _build_dependencies(self) -> str:
        """
        Builds the dependencies.
        :return: Command to load dependencies
        """
        return '' if len(self._dependencies) == 0 else 'module load {}; '.format(' '.join(self._dependencies))

    def _build_options(self, excluded_parameters: list[str] = None, delimiter: str = ' ') -> list[str]:
        """
        Builds the options string.
        :parameter delimiter: Delimiter between option and value
        :return: Options string
        """
        options = []
        for name, parameter in sorted(self._parameters.items(), key=lambda x: x[1].p_index):
            if (excluded_parameters is not None) and (name in excluded_parameters):
                continue
            if parameter.value is not None:
                options.append(parameter.option + delimiter + str(parameter.value))
            else:
                options.append(parameter.option)
        return options

    def _execute_command(self, folder: Path = None) -> None:
        """
        Executes the command.
        :return: None
        """
        if folder is None:
            folder = self._folder
        if self._command.command is None:
            raise ValueError("Command is 'None'.")
        self._command.command = self._build_dependencies() + self._command.command
        self._informs['_command'] = self._command.command
        self._command.run(folder)
        self._check_command_output()

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self.stderr != '':
            raise ToolExecutionError(f"Command execution failed (stderr: {self.stderr}).")
        elif self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    @abc.abstractmethod
    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")

    def _check_parameters(self) -> None:
        """
        Checks if the tool parameters are valid.
        :return: None
        """
        mandatory_parameters = self._tool_service.get_names_mandatory_parameter()
        for mandatory_parameter in mandatory_parameters:
            if mandatory_parameter not in self._parameters:
                raise ValueError(f"Mandatory parameter {mandatory_parameter} not set")

    def _check_input(self) -> None:
        """
        Checks if the tool input is valid.
        :return: None
        """
        for input_key, input_list in self._tool_inputs.items():
            for tool_input in input_list:
                if not isinstance(tool_input, ToolIO):
                    raise InvalidToolInputError(f"Tool input '{tool_input}' is not a ToolIO object")
                if tool_input is None:
                    raise InvalidToolInputError(f"Tool input with key {input_key} is None")
                if not tool_input.is_valid():
                    raise InvalidToolInputError(f"Invalid tool input with key {input_key}: {tool_input}")

    def _check_output(self) -> None:
        """
        Checks if the output is valid.
        :return: None
        """
        for output_key, output_list in self._tool_outputs.items():
            for tool_output in output_list:
                if tool_output is None:
                    raise ToolExecutionError(f"Tool output with key {output_key} is None")
                if not isinstance(tool_output, ToolIO):
                    raise ToolExecutionError(f"'{tool_output} {type(tool_output)}' is not a tool output object")
                if not tool_output.is_valid():
                    raise ToolExecutionError(f"Invalid tool output with key {output_key}: {tool_output}")
