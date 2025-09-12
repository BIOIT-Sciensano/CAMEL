from pathlib import Path

import pandas as pd

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class CheckM(Tool):
    """
    CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or
    metagenomes. It provides robust estimates of genome completeness and contamination by using collocated sets of genes
    that are ubiquitous and single-copy within a phylogenetic lineage. Assessment of genome quality can also be examined
    using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight sequences outside the
    expected distributions of a typical genome. CheckM also provides tools for identifying genome bins that are likely
    candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity
    within a reference genome tree.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('CheckM', '1.2.2')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_input = self.__symlink_input_files()
        dir_out = Path(self._folder) / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            'lineage_wf',
            str(dir_input),
            str(dir_out),
            '-x fasta',
        ] + self._build_options())
        self._execute_command()
        self.__set_output(self._command.stdout)

    def __symlink_input_files(self) -> Path:
        """
        Symlinks the input files.
        :return: Path to input directory
        """
        dir_input = Path(self._folder) / 'input'
        dir_input.mkdir(parents=True, exist_ok=True)
        for fasta_file in [Path(x.path) for x in self._tool_inputs['FASTA']]:
            (dir_input / f'{fasta_file.stem}.fasta').symlink_to(fasta_file)
        return dir_input

    def __set_output(self, stdout: str) -> None:
        """
        Sets the tool output for this tool.
        :param stdout: Command standard output
        :return: None
        """
        path_tsv_out = Path(self._folder) / 'output_checkm.tsv'
        output_lines = [line for line in stdout.splitlines() if not line.startswith('[') and len(line) > 0]
        with path_tsv_out.open('w') as handle:
            for line in output_lines:
                handle.write(line)
                handle.write('\n')
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]
        self._informs['results'] = pd.read_table(path_tsv_out).to_dict('records')

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the tool was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
