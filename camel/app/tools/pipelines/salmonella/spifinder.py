import json
from pathlib import Path

from camel.app.command.command import Command
from camel.app.error import InvalidToolInputError
from camel.app.error import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class SPIFinder(Tool):
    """
    BLAST-based detection tool for SPI's (Salmonella Pathogenicity Islands).
    """

    def __init__(self) -> None:
        """
        Initialize tool.
                :return: None
        """
        super().__init__('SPIFinder', '0.1')
        self._input_key = None

    def _execute_tool(self) -> None:
        """
        Executes the tool
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        db_dir = self._tool_inputs['DIR'][0].path
        self.__add_informs(db_dir)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(SPIFinder, self)._check_input()
        # check if exactly one of the three possible inputs is provided
        input_keys = [key for key in ('FASTQ', 'FASTQ_PE', 'FASTA') if key in self._tool_inputs]
        if len(input_keys) > 1:
            raise InvalidToolInputError("Too many inputs: Exactly one of FASTQ, FASTQ_PE or FASTA input is required.")
        elif len(input_keys) == 0:
            raise InvalidToolInputError("No inputs: exactly one of FASTQ, FASTQ_PE or FASTA input is required.")
        else:
            self._input_key = input_keys[0]
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DIR).")

    def __set_output(self) -> None:
        """
        Sets the name of the output files.
        The JSON output file is always data.json, therefore this needs to be hardcoded here.
        :return: None
        """
        self._tool_outputs['JSON'] = [ToolIOFile(self.folder / 'data.json')]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        if self._input_key == 'FASTQ_PE':
            inputs_str = ' '.join([
                str(self._tool_inputs[self._input_key][0].path),
                str(self._tool_inputs[self._input_key][1].path)
            ])
        else:
            inputs_str = str(self._tool_inputs[self._input_key][0].path)

        self._command.command = ' '.join([
            self._tool_command,
            f'-i {inputs_str}',
            f"-p {str(self._tool_inputs['DIR'][0].path)}"
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in command.stderr.lower() or command.exit_code != 0:
            raise ToolExecutionError(self.name, f"Error executing {self.name}: {command.stderr.strip()}")

    def __add_informs(self, db_dir: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param db_dir: Input database directory
        :return: None
        """
        self._informs['_tag'] = 'FASTQ' if 'FASTA' not in self._tool_inputs else 'FASTA'
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata file not found: {db_metadata_file}')
        with db_metadata_file.open('r') as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
