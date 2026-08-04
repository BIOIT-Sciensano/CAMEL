import json
from pathlib import Path
from typing import Any

from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.utils import fileutils

from camel.app.config import config
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class SpoTyping(Tool):
    """
    SpoTyping: fast and accurate in silico Mycobacterium spoligotyping from sequence reads or assembled contigs.

    Input:
        - FASTQ: 1 (SE) or 2 (PE) FASTQ files
        or
        - FASTA

    Output:
        - VAL_type_binary: Binary spoligotype
        - VAL_type_octal: Octal spoligotype
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('SpoTyping', None)

    @property
    def _input_key(self) -> str:
        """
        Returns the input key.
        :return: Input key
        """
        return 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command(f'{self._get_tool_command()} --version')
        self._execute_command(command, is_version_cmd=True)
        return command.stdout.splitlines()[-1].split(' ')[0]

    def _get_tool_command(self) -> str:
        """
        Returns the tool command.
        :return: String representing the command
        """
        if config.dependency_service == 'lmod':
            return self._tool_command.replace('.py', '')
        return self._tool_command

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('FASTA', 'FASTQ')):
            raise InvalidToolInputError('FASTA/Q input is required')
        for key, value in self._tool_inputs.items():
            if key == 'FASTQ' and not (0 < len(value) <= 2):
                raise InvalidToolInputError(
                    "Only 1 (SE) or 2 (PE) FASTQ inputs are supported"
                )
            if key == 'FASTA' and len(value) != 1:
                raise InvalidToolInputError("Only 1 FASTA input is supported")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Run command
        input_key = self._input_key
        command_parts = [self._get_tool_command(), *self._build_options()]
        if input_key == 'FASTQ':
            fq_in = []
            for io in self._tool_inputs['FASTQ']:
                if fileutils.is_gzipped(io.path):
                    path_out = self.folder / io.path.name.replace('.gz', '')
                    fileutils.gzip_extract(io.path, output_gz_file=path_out)
                    fq_in.append(ToolIOFile(path_out))
                else:
                    fq_in.append(io)
            command_parts.extend([str(io.path) for io in fq_in])
        else:
            io = self._tool_inputs['FASTA'][0]
            path_link = self.folder / io.path.name
            if path_link.is_symlink():
                path_link.unlink()
            path_link.symlink_to(io.path)
            command_parts.append(str(path_link))
        self._command.command = ' '.join(command_parts)
        self._execute_command()

        # Parse output
        type_binary, type_octal = self._parse_output_file(
            self.folder / fileutils.make_valid(self.get_param_value('output_basename'))
        )

        # Set output
        self._tool_outputs['VAL_type_binary'] = [ToolIOValue(type_binary)]
        self._tool_outputs['VAL_type_octal'] = [ToolIOValue(type_octal)]
        self._tool_outputs['LOG'] = [
            ToolIOFile(
                Path(self._folder) / f'{self.get_param_value("output_basename")}.log'
            )
        ]
        self._informs['metadata'] = self._extract_metadata(type_octal)

    def _parse_output_file(self, output_file: Path) -> tuple[str, str]:
        """
        Parses the output file.
        :param output_file: Path to output file
        :return: Spoligotype (Binary), Spoligotype (Octal)
        """
        if not output_file.exists():
            raise ToolExecutionError(self.name, "Output file not found")
        with output_file.open('r') as handle:
            try:
                _, type_binary, type_octal = handle.readlines()[-1].strip().split('\t')
                logger.info(f'Detected spoligotype (octal): {type_octal}')
                return type_binary, type_octal
            except (IndexError, ValueError):
                raise ToolExecutionError(self.name, "Output file has an invalid format")

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output to check if the tool executed successfully.
        :return: None
        """
        if command.exit_code != 0:
            last_line = command.stderr.splitlines()[-1]
            if last_line.startswith('urllib2.URLError'):
                logger.warning('Could not contact SITVIT server')
            else:
                raise ToolExecutionError(self.name, last_line)

    def _extract_metadata(self, type_octal: str) -> dict[str, Any]:
        """
        Extracts the metadata for the detected Spoligotype.
        :return: Spoligotype metadata
        """
        # Get location of metadata file
        command = Command(f'{self._build_dependencies()} echo $SPOTYPING_METADATA')
        command.run(self.folder)
        metadata_path = command.stdout.strip()
        if not Path(metadata_path).exists():
            raise FileNotFoundError("No spoligotype metadata found")

        # Parse info
        with open(metadata_path) as handle:
            metadata = json.load(handle)

        # Extract metadata
        keys = ('SIT', 'geo', 'label', 'total')
        if type_octal in metadata:
            return {k: metadata[type_octal][k] for k in keys}
        else:
            return {k: 'NA' for k in keys}
