import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Union, Optional, Any

from snakemake.io import Wildcards

from camel.app import loggers
from camel.app.components import toolutils
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


@dataclass(frozen=True)
class StepOutput:
    """
    This class is used to keep an output of a step.
    """
    rule_name: str
    type: str
    key: str
    index: int
    hash: str
    wildcards: Optional[Wildcards] = None

    @property
    def wildcards_json(self) -> Union[None, str]:
        """
        Returns the wildcards as a JSON string.
        :return: Wildcards as string
        """
        if self.wildcards is None:
            return None
        return json.dumps({k: v for k, v in self.wildcards.items()})


class Step:
    """
    This class represents a step in a Snakemake pipeline. It executes a single tool.
    """

    def __init__(self, rule_name: str, tool: Tool, wildcards: Wildcards | Any = None, dir_: Optional[Path]= None) -> None:
        """
        Initializes a step.
        :param rule_name: Name of the rule in Snakemake.
        :param tool: Tool object to execute.
        :param dir_: Working directory.
        :param wildcards: Wildcards object from snakemake.
        :return: None
        """
        self._name = rule_name
        self._tool = tool
        # self._step_inputs = {}
        # self._input_informs = {}
        self._dir = Path(dir_).absolute()
        self._wildcards = wildcards

    @property
    def name(self) -> str:
        """
        Returns the name of this step.
        :return: Name
        """
        return self._name

    @property
    def informs(self) -> dict:
        """
        Returns the informs of this step.
        :return: Informs
        """
        return self._tool.informs

    # def add_inputs(self, dict_: dict) -> None:
    #     """
    #     Adds the inputs to the step
    #     :param dict_: Dictionary with input objects
    #     :return: None
    #     """
    #     self._step_inputs = dict_

    # def add_informs(self, dict_: dict) -> None:
    #     """
    #     Adds informs to the step.
    #     :param dict_: Dictionary with the informs
    #     :return: None
    #     """
    #     self._input_informs = dict_
    #     logger.info(f"Inform added: {dict_}")

    def run(self) -> None:
        """
        Runs this step.
        :return: None
        """
        loggers.attach_step_handler(self._dir / 'logs', logging.DEBUG)
        with Path(self._dir / 'logs' / 'rulename.txt').open('a') as handle:
            handle.write(self.name)
            handle.write('\n')
        logger.info(f'Running step: {self.name}')
        # self._tool.add_input_files(self._step_inputs)
        # self._tool.add_input_informs(self._input_informs)
        logger.debug(f"Tool inputs: {self._tool.tool_inputs}")
        logger.info(f"Default parameters loaded: {toolutils.show_parameters(self._tool)}")
        dir_tool = self._dir / 'work'
        dir_tool.mkdir(exist_ok=True)
        self._tool.run(dir_tool)
        logger.info(f'Step output: {list(self._tool.tool_outputs.items())}')
        logger.info(f'Step informs: {list(self._tool.informs.items())}')
        self._log_outputs()
        loggers.detach_step_handlers()

    def run_step(self) -> None:
        """
        Runs the current step.
        :return: None
        """
        logger.warning('.run_step() is deprecated, please use .run() instead.')
        self.run()

    def _log_outputs(self) -> None:
        """
        Logs the outputs in the database.
        :return: None
        """
        logger.info(f"Logging output for step '{self.name}'")
        for key, io_list in self._tool.tool_outputs.items():
            for index, io_out in enumerate(io_list):
                if not io_out.is_logged:
                    continue
                step_output = StepOutput(self._name, io_out.type_name, key, index, io_out.hash, self._wildcards)
                logger.debug(f'Output log: {step_output}')
