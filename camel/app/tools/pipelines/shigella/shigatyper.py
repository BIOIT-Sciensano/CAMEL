from pathlib import Path

import pandas as pd
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError, ToolExecutionError
from camel.app.core.tool import Tool
from camel.app.loggers import logger


class ShigaTyper(Tool):
    """
    ShigaTyper does Shigella/EIEC identification and serotyping based on WGS data.
    """

    def __init__(self) -> None:
        """
        Initializes the ShigaTyper tool.
        """
        super().__init__('ShigaTyper', '2.0.5')

    def _check_input(self) -> None:
        """
        Checks whether the provided input is valid:
        - Illumina paired-end reads are the only required input
        :param: fastq_pe : List of forward and reverse fastq files
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('FASTQ_PE', 'FASTQ_SE')):
            raise InvalidToolInputError('FASTQ_PE or FASTQ_SE input is required.')
        super()._check_input()

    def __build_command(self, sample_name: str) -> None:
        """
        Concatenates required parameters and options to build the command.
        :param sample_name: Sample ID
        :return: None
        """
        # Construct the read arguments
        if 'FASTQ_PE' in self._tool_inputs:
            read_args = f"--R1 {self._tool_inputs['FASTQ_PE'][0].path} --R2 {self._tool_inputs['FASTQ_PE'][1].path}"
        else:
            read_args = f"--SE {self._tool_inputs['FASTQ_SE'][0].path}"

        # Create the full command
        self._command.command = ' '.join([
            self._tool_command,
            read_args,
            f'--name {sample_name}',
            *self._build_options()
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Run the command
        self.__build_command('shigatyper_out')
        self._execute_command()

        # Create dummy output for isolates not covered by the tool
        create_dummy_output = False
        if 'Checkpoint 1 failed' in self._command.stderr:
            create_dummy_output = True

        # Collect the output
        for key, basename in zip(('TSV', 'TSV_HITS'), ('shigatyper_out.tsv', 'shigatyper_out-hits.tsv')):
            path_out = self.folder / basename
            if not path_out.exists():
                if key == 'TSV_HITS' and create_dummy_output is True:
                    logger.info(f'Creating dummy output file: {path_out}')
                    path_out.touch()
                else:
                    raise ToolExecutionError(self.name, f'{path_out} not generated ({key})')
            self._tool_outputs[key] = [ToolIOFile(path_out)]

        # Parse TSV output file
        self._parse_tsv(self._tool_outputs['TSV'][0].path)

    def _parse_tsv(self, path_tsv: Path) -> None:
        """
        Parses the output TSV file and stores the results in the informs.
        :param path_tsv: Path to output file
        :return: None
        """
        data_serotype = pd.read_table(path_tsv)
        output_dict = data_serotype.to_dict('records')[0]
        self._informs['species'] = output_dict['prediction']
