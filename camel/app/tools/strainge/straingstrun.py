from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class StrainGSTRun(Tool):
    """
    StrainGST (Strain Genome Search tool) is a tool to find close reference genomes for strains present in a sample.
    StrainGST run is used to query a FASTQ file against a database. Both should already be formatted in HDF5 format.
    """

    def __init__(self) -> None:
        """
        Initializes straingst run.
                :return: None
        """
        super().__init__('StrainGST run', '1.3.9')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'DB_HDF5' not in self._tool_inputs:
            raise InvalidToolInputError('HDF5 database is required')
        if 'HDF5' not in self._tool_inputs:
            raise InvalidToolInputError('input HDF5 file is required')
        super()._check_input()

    def _build_command(self, hdf5_input: Path, db_input: Path) -> None:
        """
        Builds the command to run StrainGST run.
        :param hdf5_input: Path to the input file in HDF5 format
        :param db_input: Path to the database in HDF5 format
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            f'{db_input}',
            f'{hdf5_input}',
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def _set_output(self) -> None:
        """
        Collects the tool output.
        :return: None
        """
        prefix_string = self._parameters["output_prefix"].value
        self._tool_outputs['TSV_STATS'] = \
            [ToolIOFile(self.folder / f'{prefix_string}.stats.tsv')]
        self._tool_outputs['TSV_STRAINS'] = \
            [ToolIOFile(self.folder / f'{prefix_string}.strains.tsv')]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        hdf5_input = Path(str(self._tool_inputs['HDF5'][0]))
        db_input = Path(str(self._tool_inputs['DB_HDF5'][0]))
        self._build_command(hdf5_input, db_input)
        self._execute_command()
        self._set_output()
