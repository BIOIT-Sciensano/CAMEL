from pathlib import Path

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class BAMAddCustomTag(Tool):
    """
    Tool that adds a custom tag to a BAM file.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('BAM add custom tag', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidToolInputError('BAM input is required')
        super()._check_input()

    def __extract_header(self, bam: Path) -> str:
        """
        Extracts the original header.
        :param bam: Input BAM file
        :return: Header content
        """
        command = Command(f'{self._build_dependencies()} samtools view -H {bam}')
        command.run(self.folder)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error extracting header from BAM file: {bam}')
        return command.stdout

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Add custom tag and save header
        header = self.__extract_header(self._tool_inputs['BAM'][0].path)
        header += f"@CO\t{self._parameters['name'].value}:{self._parameters['value'].value}\n"
        path_header_updated = self.folder / 'header_updated.txt'
        with path_header_updated.open('w') as handle:
            handle.write(header)
        logger.info(f'Updated header file created:{path_header_updated}')

        # Create output BAM file with custom tag
        path_out = self.folder / self._parameters['output'].value
        self._command.command = ' '.join([
            'samtools reheader',
            str(path_header_updated),
            str(self._tool_inputs['BAM'][0].path),
            f'> {path_out}'
        ])
        self._execute_command()

        # Set the output
        self._tool_outputs['BAM'] = [ToolIOFile(path_out)]
