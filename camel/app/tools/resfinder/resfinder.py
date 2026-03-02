import json
from pathlib import Path

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class ResFinder(Tool):
    """
    ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.
    """

    def __init__(self) -> None:
        """
        Initializes this tool
        :return: None
        """
        super().__init__('ResFinder', '4.7.2')

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if any(key in self._tool_inputs for key in ('FASTA' or 'FASTQ_PE' or 'FASTQ_SE')):
            raise InvalidToolInputError('FASTA or FASTQ input is required')
        if not ('acquired' or 'point' in self._tool_inputs):
            raise InvalidToolInputError('Either "acquired" or "point" is required')
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database input is required (DIR)")
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run resfinder.
        :return: None
        """
        if 'FASTA' in self._tool_inputs:
            input_str = f'--inputfasta {self._tool_inputs["FASTA"][0].path}'
        elif 'FASTQ_SE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ"][0].path}'
        elif 'FASTQ_PE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ_PE"][0].path} ' \
                        f'{str(self._tool_inputs["FASTQ_PE"][1].path)}'
        else:
            raise ValueError('Invalid tool input')
        self._command.command = ' '.join([
            self._tool_command, input_str, '-db_point', str(self._tool_inputs['DIR'][0].path / 'pointfinder'),
            '-db_res', str(self._tool_inputs['DIR'][0].path / 'resfinder'), *self._build_options()])

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
        """
        dir_out = self.folder / self._parameters['output_path'].value
        self._tool_outputs['TSV_pheno_general'] = [ToolIOFile(dir_out / Path('pheno_table.txt'))]
        try:
            self._informs['species'] = self._parameters['species'].value
        except KeyError:
            self._informs['species'] = ''

        if 'acquired' in self._parameters:
            self._tool_outputs['TSV_genes'] = [ToolIOFile(dir_out / Path('ResFinder_results_tab.txt'))]
        if 'point' in self._parameters:
            self._tool_outputs['TSV_point'] = [ToolIOFile(dir_out / Path('PointFinder_results.txt'))]
            try:
                self._tool_outputs['TSV_pheno_species'] = [ToolIOFile(next(dir_out.glob('pheno_table_*.txt')))]
            except StopIteration:
                self._tool_outputs['TSV_pheno_species'] = []

    def __collect_db_version(self) -> None:
        """
        Collect the db version from the dbupdate json output file.
        :return: None
        """
        with open(self._tool_inputs['DIR'][0].path / 'resfinder' / 'db_update_info.json') as handle:
            db_version_json_resfinder = json.load(handle)
        with open(self._tool_inputs['DIR'][0].path / 'pointfinder' / 'db_update_info.json') as handle:
            db_version_json_pointfinder = json.load(handle)
        self._informs['db_version_resfinder'] = db_version_json_resfinder['last_update_date']
        self._informs['db_version_pointfinder'] = db_version_json_pointfinder['last_update_date']
        path_version_file = Path(self._tool_inputs['DIR'][0].path, 'resfinder', 'VERSION')
        if path_version_file.exists():
            with path_version_file.open() as handle:
                db_version = handle.readline().strip()
        self._informs['db_version_name'] = db_version if path_version_file.exists() else 'n/a'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
        self._set_informs()
        self.__collect_db_version()

    def _set_informs(self) -> None:
        """
        Collects the tool informs.
        :return: None
        """
        self._informs['parameters'] = {
            'threshold': float(self._parameters['threshold'].value) if 'threshold' in self._parameters else 'default',
            'min_cov': float(self._parameters['min_cov'].value) if 'min_cov' in self._parameters else 'default',
        }
