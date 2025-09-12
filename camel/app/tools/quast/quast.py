from pathlib import Path

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.error import InvalidToolInputError, ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class Quast(Tool):
    """
    QUAST evaluates genome assemblies. QUAST works both with and without a reference genome. The tool accepts multiple
    assemblies, thus is suitable for comparison.
    """

    def __init__(self) -> None:
        """
        Initialize tool
        :return: None
        """
        super().__init__('quast', '5.2.0')

    def _execute_tool(self) -> None:
        """
        Runs QUAST.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA is required
        - FASTA_Ref, TSV_Gene, and TSV_Operon are optional
        - Only one input file allowed for FASTA_Ref, TSV_Gene, and TSV_Operon, multiple files allowed for FASTA
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError(
                f'QUAST required FASTA input is missing: {self._tool_inputs!r}')
        for key, values in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTA_Ref', 'TSV_Gene', 'TSV_Operon', 'BAM', 'BAM_Ref', 'FASTQ_PE',
                           'FASTQ_nanopore', 'GFF3_Ref']:
                raise InvalidToolInputError(
                    f'Illegal input key given for QUAST: {self._tool_inputs!r}')
            if key in ['FASTA_Ref', 'TSV_Gene', 'TSV_Operon'] and len(values) > 1:
                raise InvalidToolInputError(
                    f'Too many input files given for QUAST: {self._tool_inputs!r}'
                )
        super()._check_input()

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command.
        :return: None
        """
        options_string = ' '.join(self._build_options() + [f'-o {self._folder}', '--debug'])
        input_string = self.__build_input_string()
        self._command.command = ' '.join([self._tool_command, options_string, input_string])

    def __build_input_string(self) -> str:
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        inputs = []
        if 'FASTA_Ref' in self._tool_inputs:
            inputs.append(f"-r {self._tool_inputs['FASTA_Ref'][0].path}")
        if 'GFF3_Ref' in self._tool_inputs:
            inputs.append(f"--features {self._tool_inputs['GFF3_Ref'][0].path}")
        if 'TSV_Gene' in self._tool_inputs:
            inputs.append(f"-G {self._tool_inputs['TSV_Gene'][0].path}")
        if 'TSV_Operon' in self._tool_inputs:
            inputs.append(f"-O {self._tool_inputs['TSV_Operon'][0].path}")
        if 'BAM' in self._tool_inputs:
            inputs.append(f"--bam {self._tool_inputs['BAM'][0].path}")
        if 'BAM_Ref' in self._tool_inputs:
            inputs.append(f"--ref-bam {self._tool_inputs['BAM_Ref'][0].path}")
        if 'FASTQ_PE' in self._tool_inputs:
            inputs.extend([
                f"--pe1 {self._tool_inputs['FASTQ_PE'][0].path}", f"--pe2 {self._tool_inputs['FASTQ_PE'][1].path}"])
        if 'FASTQ_nanopore' in self._tool_inputs:
            inputs.append(f"--nanopore {self._tool_inputs['FASTQ_nanopore'][0].path}")
        for item in self._tool_inputs['FASTA']:
            inputs.append(str(item.path))
        return ' '.join(inputs)

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        for line in command.stderr.splitlines():
            if 'ERROR' in line:
                if 'ERRORs: 0' not in line:
                    raise ToolExecutionError(self.name, f"Command execution failed (stderr: {command.stderr}).")
        toolutils.check_tool_execution(self, command, exit_code=0)

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        output_keys = ['HTML', 'TEX', 'TSV', 'TXT']
        for key in output_keys:
            self._tool_outputs[key] = [ToolIOFile(Path(f'{self.folder / "report"}.{key.lower()}'))]

        # for icarus browser
        icarus_output_keys = {
            'HTML_icarus': 'icarus.html',
            'HTML_alignment_viewer': 'icarus_viewers/alignment_viewer.html',
            'HTML_contig_size_viewer': 'icarus_viewers/contig_size_viewer.html'
        }
        for key, value in icarus_output_keys.items():
            if key == 'HTML_alignment_viewer' and 'FASTA_Ref' not in self._tool_inputs:
                # skip HTML_alignment_viewer when no reference genome is provided
                continue
            path_out = self.folder / value
            if not path_out.exists():
                logger.warning(f'HTML output not found: {path_out}')
                continue
            self._tool_outputs[key] = [ToolIOFile(self.folder / value)]

        # Optional analyses
        if 'glimmer' in self._parameters:
            self._tool_outputs['GFF'] = [ToolIOFile(next((self.folder / 'predicted_genes').glob('*.gff')))]
        if 'conserved_genes_finding' in self._parameters:
            self._tool_outputs['TXT_busco'] = [ToolIOFile(next((self.folder / 'busco_stats').glob(
                'short_summary_*.txt')))]

        # Reference genome
        if 'FASTA_Ref' in self._tool_inputs:
            self._informs['ref'] = self._tool_inputs['FASTA_Ref'][0].path.name
