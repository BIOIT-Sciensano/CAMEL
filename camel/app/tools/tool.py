import abc
import inspect
from pathlib import Path
from typing import Optional, Union, Any
from collections import OrderedDict

import pydantic
import yaml

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.error import InvalidParameterError
from camel.app.error import ToolExecutionError
from camel.app.io.toolio import ToolIO
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.parameter.parameter import Parameter


class Tool(metaclass=abc.ABCMeta):
    """
    Base class for tool classes.
    """

    def __init__(self, name: str, version: str | None, camel: Optional[Camel] = None) -> None:
        """
        Initializes a tool.
        :param name: Tool name
        :param version: Tool version
        :return: None
        """
        self._name: str = name

        # Read YAML data
        tool_data = self._read_tool_yml()
        self._tool_command = tool_data['tool_command']
        self._dependencies = tool_data['dependencies']

        # Command & directory
        self._command = Command()
        self._folder: Optional[Path] = None

        # Get the tool version
        self._version: str = version if version is not None else self.get_version()

        # Setup inputs / outputs
        logger.debug(f"Initializing tool: {self._name} {self._version}")
        self._tool_inputs: dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]] = {}
        self._tool_outputs: dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]] = {}
        self._informs: dict[str, Any] = {'_name': self.name, '_version': self._version}
        self._input_informs = {}

        # Parameters
        self._param_data: dict[str, Any] = tool_data['parameters'] # All available parameters
        self._params: OrderedDict[str, Parameter] = self.get_default_params() # Currently active parameters

        # To deprecate
        self._camel = camel

    def get_version(self) -> str:
        """
        This method can be implemented by subclasses to dynamically retrieve the tool version.
        :return: Tool version (as a string)
        """
        if self._version is None:
            raise RuntimeError('get_version should be implemented when no version is provided')
        return self._version

    @property
    def _parameters(self) -> dict:
        """
        Placeholder for accessing parameters (for legacy reasons).
        :return: Parameter dict
        """
        return self._params

    def param_data(self) -> dict:
        """
        Returns the parameter data.
        :return: Parameter data
        """
        return self._param_data

    @property
    def name(self) -> str:
        """
        Returns the name of this tool.
        :return: Name
        """
        return self._name

    @property
    def name_full(self) -> str:
        """
        Returns the full name of this tool (including tool version).
        :return: Name
        """
        return f'{self.name} {self.version}'

    @property
    def version(self) -> str:
        """
        Returns the tool version.
        :return: Version
        """
        return self._version

    @property
    def tool_outputs(self) -> dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]]:
        """
        Returns the tool outputs.
        :return: Tool outputs
        """
        return self._tool_outputs

    @property
    def tool_inputs(self) -> dict[str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]]:
        """
        Returns the tool input.
        :return: Tool inputs
        """
        return self._tool_inputs

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
    def params(self) -> dict[str, Parameter]:
        """
        Returns the dictionary with currently active parameters.
        :return: Parameter dictionary
        """
        return self._params

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
        for param_key, new_value in kwargs.items():
            parameter = self.get_param(param_key)
            if not parameter:
                raise InvalidParameterError(f"{self._name} has no parameter '{param_key}'")

            # Value is False -> remove the parameter
            if new_value is False:
                if param_key not in self._params:
                    logger.warning(f"Cannot disable parameter '{param_key}' (not present in parameters)")
                    continue
                logger.info(f"Disabling parameter: {param_key}")
                self._params.pop(param_key)
                continue

            # Flag -> add it to the parameter dict
            if parameter.flag and new_value is True:
                self._params[param_key] = parameter
                continue

            # Value change existing
            if param_key in self.params:
                old_value = self._params[param_key].value
                logger.debug(f"Parameter '{param_key}' value '{old_value}' changed to '{new_value}'")
                self._params[param_key].value = new_value
                continue

            # Set parameter
            parameter.value = new_value
            self._params[param_key] = parameter

    def clear_parameters(self) -> None:
        """
        Clears all the parameters of the given tool.
        :return: None
        """
        logger.info(f"Removing {len(self._params)} parameters")
        self._params.clear()

    def run(self, folder: Path = Path.cwd()) -> None:
        """
        Runs this tool.
        :param folder: Folder to run the tool in.
        :return: None
        """
        self._folder = folder
        logger.info(f'Running tool {self.name}')
        logger.info(f'Working directory: {self._folder}')
        logger.info(f'Tool parameters: {toolutils.show_parameters(self)}')
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
        for name, parameter in sorted(self._params.items(), key=lambda x: x[1].p_index):
            if (excluded_parameters is not None) and (name in excluded_parameters):
                continue
            if parameter.value is True or parameter.value is None:
                if parameter.flag is False:
                    logger.warning(f'Consider changing {name} into a flag parameter')
                options.append(parameter.option)
            else:
                options.append(parameter.option + delimiter + str(parameter.value))
        return options

    def _execute_command(self, command: Command | None = None, dir_: Path | None = None, is_version_cmd: bool = False) -> None:
        """
        Executes the given command.
        :param command: Command to execute. Defaults to self._command.
        :param dir_: Directory in which to execute the command.
        :param is_version_cmd: Boolean to indicate is this is a version command
        :return: None
        """
        # Setup command
        if command is None:
            command = self._command
        if command.command is None:
            raise ValueError("Command should not be None")

        # Setup working directory
        if dir_ is None:
            dir_ = self._folder if not is_version_cmd else Path.cwd()

        # Execute the command
        if not is_version_cmd:
            self._informs['_command'] = self._build_dependencies() + command.command
        command.run(dir_, prefix=self._build_dependencies())
        self._check_command_output(command)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        toolutils.check_tool_execution(self, command)

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
        mandatory = [key for key, value in self._param_data.items() if value.get('mandatory', False) is True]
        for key in mandatory:
            if key in self._params:
                continue
            raise InvalidParameterError(f"Mandatory parameter '{key}' not set")

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
                    raise ToolExecutionError(self.name, f"Tool output with key {output_key} is None")
                if not isinstance(tool_output, ToolIO):
                    raise ToolExecutionError(self.name, f"'{tool_output} {type(tool_output)}' is not a ToolIO object")
                if not tool_output.is_valid():
                    raise ToolExecutionError(self.name, f"Invalid tool output with key {output_key}: {tool_output}")

    def _get_tool_data_yml_path(self) -> Path:
        """
        Returns the path to the tool YAML file.
        :return: Path
        """
        path_yaml = Path(inspect.getfile(self.__class__).replace('.py', '.yml'))
        if not path_yaml.is_file():
            raise FileNotFoundError(f"Tool data file for '{self.name}' not found: {path_yaml}")
        return path_yaml

    def _read_tool_yml(self) -> dict[str, Any]:
        """
        Reads the tool YAML file.
        :return: Tool data as a dictionary
        """
        path_yml = self._get_tool_data_yml_path()
        with path_yml.open() as handle:
            return yaml.safe_load(handle)

    def get_default_params(self) -> OrderedDict[str, Parameter]:
        """
        Returns the default parameters for this tool.
        :return: Default parameters
        """
        param_dict = OrderedDict()
        for p_name, p_data in sorted(self._param_data.items(), key=lambda x: x[1].get('p_index', 0)):
            if p_data.get('default', False) is False:
                continue
            try:
                param_dict[p_name] = Parameter(**{'name': p_name, **p_data})
            except pydantic.ValidationError as err:
                logger.error(f'Invalid parameter data for {p_name}: {err}')
                raise err
        return param_dict

    def get_param(self, key: str) -> Optional[Parameter]:
        """
        Returns the parameter with the given key.
        :param key: Parameter key
        :return: Parameter (if available, None otherwise)
        """
        try:
            return Parameter(**{'name': key, **self._param_data[key]})
        except KeyError:
            return None
        except pydantic.ValidationError as err:
            logger.error(f'Invalid parameter data for {key}: {err}')
            raise err

    def get_param_value(self, key: str) -> Any:
        """
        Returns the value for the parameter with the given key.
        :param key: Parameter key
        :return: Parameter value
        """
        if key not in self._param_data:
            raise InvalidParameterError(f'Tool {self.name} has no parameter {key}')
        param = self.get_param(key)
        if param.flag:
            return key in self._params
        if key in self._params:
            return self._params[key].value
        return None
