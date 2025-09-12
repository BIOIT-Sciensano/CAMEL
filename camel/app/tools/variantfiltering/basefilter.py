from pathlib import Path

import abc

from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class BaseFilter(Tool, metaclass=abc.ABCMeta):
    """
    Base class for variant filters.
    """

    def __init__(self, *args) -> None:
        """
        Initializes this filter.
        :param args: Initialization arguments
        :return: None
        """
        super().__init__(*args)
        self._output_file = None

    @property
    @abc.abstractmethod
    def full_name(self) -> str:
        """
        The full name for the filter.
        :return: Full name
        """
        pass

    @property
    @abc.abstractmethod
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        raise NotImplementedError()

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'VCF_GZ' not in self._tool_inputs:
            raise InvalidToolInputError("No compressed VCF input found")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        nb_of_variants_pre = VCFUtils.count_variants(self._tool_inputs['VCF_GZ'][0].path)
        self._apply_filter()
        nb_of_variants_post = VCFUtils.count_variants(self.output_path)
        logger.info(f'{nb_of_variants_post}/{nb_of_variants_pre} variants passed {self._name} filtering')
        self._informs['variants_in'] = nb_of_variants_pre
        self._informs['variants_out'] = nb_of_variants_post
        self._informs['full_name'] = self.full_name
        self._informs['description'] = self.description
        self._tool_outputs['VCF_GZ'] = [ToolIOFile(self.output_path)]

    @abc.abstractmethod
    def _apply_filter(self) -> None:
        """
        Applies the filter to the input file.
        :return: None
        """
        pass

    @property
    def output_path(self) -> Path:
        """
        Returns the path to the output file.
        :return: Path
        """
        return Path(self._folder) / self._parameters['output_filename'].value

    def _get_soft_filter_options(self) -> list[str]:
        """
        Returns the parts of the command related to the soft filtering.
        If soft filtering is not enabled an empty list is returned.
        :return: List of options
        """
        if 'soft_filter' not in self._parameters:
            return []
        filter_name = self.full_name.replace(' ', '_').lower()
        return [f"{self._parameters['soft_filter'].option} '{filter_name}'"]
