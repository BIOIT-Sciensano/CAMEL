from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolsbase import BcftoolsBase


class BcftoolsFilter(BcftoolsBase):
    """
    Filtering of VCF/BCF files using fixed thresholds.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('bcftools filter', '1.17')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('VCF', 'VCF_GZ')):
            raise InvalidToolInputError("VCF/VCF_GZ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._tool_outputs[self._get_output_key()] = [ToolIOFile(self._get_output_path())]

    def _build_command(self) -> None:
        """
        Builds the command that is called.
        :return: None
        """
        command_parts = [self._tool_command]
        if 'BED_include' in self._tool_inputs:
            command_parts.append(f"--targets-file {self._tool_inputs['BED_include'][0].path}")
        elif 'BED_exclude' in self._tool_inputs:
            command_parts.append(f"--targets-file ^{self._tool_inputs['BED_exclude'][0].path}")
        command_parts.extend(self._build_options(excluded_parameters=['invert_targets']))
        command_parts.append(str(next(
            self._tool_inputs[k][0].path for k in ('VCF', 'VCF_GZ') if k in self._tool_inputs.keys())))
        self._command.command = ' '.join(command_parts)
