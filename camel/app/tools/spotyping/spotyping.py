import json
from pathlib import Path
from typing import Any

from camel.app.command.command import Command
from camel.app.error import InvalidToolInputError, ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


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
        super().__init__('SpoTyping', '2.1')
        self._input_key = None

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('FASTA', 'FASTQ')):
            raise InvalidToolInputError('FASTA/Q input is required')
        for key, value in self._tool_inputs.items():
            if key == 'FASTQ' and not (0 < len(value) <= 2):
                raise InvalidToolInputError("Only 1 (SE) or 2 (PE) FASTQ inputs are supported")
            if key == 'FASTA' and len(value) != 1:
                raise InvalidToolInputError("Only 1 FASTA input is supported")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Run command
        self._input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        path_inputs = self._symlink_input()
        self._command.command = ' '.join([
            self._tool_command,
            *[str(f) for f in path_inputs],
            * self._build_options()
        ])
        self._execute_command()

        # Parse output
        type_binary, type_octal = SpoTyping._parse_output_file(
            Path(self._folder) / self._parameters['output_basename'].value)

        # Set output
        self._tool_outputs['VAL_type_binary'] = [ToolIOValue(type_binary)]
        self._tool_outputs['VAL_type_octal'] = [ToolIOValue(type_octal)]
        self._tool_outputs['LOG'] = [ToolIOFile(Path(self._folder) / f'{self.get_param_value("output_basename")}.log')]
        self._informs['metadata'] = self._extract_metadata(type_octal)

    @staticmethod
    def _parse_output_file(output_file: Path) -> tuple[str, str]:
        """
        Parses the output file.
        :return: Spoligotype (Binary), Spoligotype (Octal)
        """
        if not output_file.exists():
            raise ToolExecutionError('SpoTyping', "Output file not found")
        with output_file.open('r') as handle:
            try:
                _, type_binary, type_octal = handle.readlines()[-1].strip().split('\t')
                return type_binary, type_octal
            except IndexError:
                raise ToolExecutionError('SpoTyping', "Output file has an invalid format")

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output to check if the tool executed successfully.
        :return: None
        """
        if self._command.exit_code != 0:
            last_line = self._command.stderr.splitlines()[-1]
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
        command.run(self._folder)
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

    def _symlink_input(self) -> list[Path]:
        """
        Symlinks the input file(s).
        :return: Paths to symlinked input files
        """
        path_inputs = [self._folder / f.path.name for f in self._tool_inputs[self._input_key]]
        for path_new, path_old in zip(path_inputs, self._tool_inputs[self._input_key]):
            if path_new.is_symlink():
                path_new.unlink()
            path_new.symlink_to(path_old.path)
        return path_inputs
