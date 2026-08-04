import json
import re
from pathlib import Path

import pandas as pd
from camelcore.app.command import Command
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class Mykrobe(Tool):
    """
    Mykrobe performs antibiotic resistance prediction and genotyping for Mycobacterium tuberculosis,
    Staphylococcus aureus, Shigella sonnei and Salmonella typhi.
    """

    def __init__(self) -> None:
        """
        Initializes Mykrobe.
        """
        super().__init__('mykrobe', '0.13.0')

    def _execute_tool(self) -> None:
        """
        Runs Mykrobe.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        input_folder = self._tool_inputs['DIR'][0].path
        self.__add_informs(input_folder)
        self._parse_csv(self._tool_outputs['CSV'][0].path)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        input_types = {'FASTQ_PE', 'FASTA', 'FASTQ_SE'}

        super()._check_input()
        if 'SPECIES' not in self._tool_inputs:
            raise InvalidToolInputError(
                "Species (i.e. 'sonnei', 'staph', 'tb' or 'typhi') needs to be specified")
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("Database path needs to be specified")
        if len(input_types.intersection(set(self._tool_inputs))) != 1:
            raise InvalidToolInputError(
                "One input type (i.e. 'FASTQ_PE', 'FASTA' or 'FASTQ_SE') is required and accepted")

    @staticmethod
    def __prepare_input_str(input_type: dict) -> str:
        """
        Builds input based on format (i.e. 'FASTQ_PE', 'FASTA' or 'FASTQ_SE').
        :return: None
        """
        if 'FASTQ_PE' in input_type:
            return f"{str(input_type['FASTQ_PE'][0].path)} {str(input_type['FASTQ_PE'][1].path)}"
        elif 'FASTA' in input_type:
            return str(input_type['FASTA'][0].path)
        elif 'FASTQ_SE' in input_type:
            return str(input_type["FASTQ_SE"][0].path) + ' --ont'
        raise ValueError(f'Invalid input type: {input_type}')

    def __build_command(self) -> None:
        """
        Builds the command line call to execute Mykrobe.
        :return: None
        """
        input_str = self.__prepare_input_str(self._tool_inputs)
        species_flag = self._tool_inputs['SPECIES'][0]
        sample_id = re.sub('_[1-2]', '', Path(input_str).stem.replace('.fastq', ''))

        self._command.command = ' '.join([
            self._tool_command, 'predict',
            f'--sample {sample_id}',
            f'--species {species_flag}',
            '--format csv',
            f'--seq {input_str}',
            *self._build_options()
        ])

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(self.folder / str(self._parameters['output_filename'].value))]

    def __add_informs(self, input_folder: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param input_folder: Input database directory
        :return: None
        """
        path_metadata = input_folder / 'db_update_info.json'
        if not path_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {path_metadata}')
        with path_metadata.open('r') as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        self._informs['db_version'] = metadata['last_update_date']

    def _parse_csv(self, path_csv: Path) -> None:
        """
        Parses the output CSV file and stores the results in the informs.
        :param path_csv: Path to output file
        :return: None
        """
        data = pd.read_csv(path_csv)
        self._informs['phylo_group'] = data['phylo_group'][0]
        self._informs['species'] = data['species'][0]
        self._informs['lineage'] = data['lineage'][0]
        self._informs['drug_susceptibility'] = data.iloc[:, 1:5].values.tolist()
