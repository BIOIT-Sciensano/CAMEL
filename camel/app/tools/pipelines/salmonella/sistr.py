import json
from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Sistr(Tool):
    """
    Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and
    cgMLST gene alleles using BLAST.
    """
    def __init__(self) -> None:
        """
        Initialize tool.
                :return: None
        """
        super().__init__('SISTR', '1.1.1')

    def _execute_tool(self):
        """
        Execute the tool.
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
        super()._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DIR).")

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['JSON'] = [ToolIOFile(self.folder / self._parameters['output_filename'].value)]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '--output-format json',
            '--use-full-cgmlst-db',
            '--qc',
            '--verbose',
            *self._build_options(),
            str(self._tool_inputs['FASTA'][0].path)
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        We don't check "if 'error' in self.stderr.lower()" here because some small warnings are wrongfully displayed by
        the tool as errors.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __add_informs(self, db_dir: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param db_dir: Input database directory
        :return: None
        """
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata not found: {db_metadata_file}')
        with db_metadata_file.open() as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        self._informs['db_path'] = str(db_dir)
