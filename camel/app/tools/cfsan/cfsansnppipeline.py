from importlib.resources import files
from pathlib import Path

import os
import shutil

from camel.app.core.command import Command
from camel.app.core.utils import toolutils, fileutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class CfsanSnpPipeline(Tool):
    """
    Runs the SNP pipeline developed by CFSAN.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('CFSAN SNP Pipeline', '2.2.1')

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("No reference FASTA input found.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidToolInputError("Only one reference FASTA file is supported.")
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidToolInputError("No FASTQ input found.")
        if 'VAL_name' not in self._tool_inputs:
            raise InvalidToolInputError("No sample name input found ('VAL_name')")
        if len(self._tool_inputs['FASTQ']) % 2 != 0:
            raise InvalidToolInputError("Only paired end input is supported (2 entries per sample)")
        if len(self._tool_inputs['VAL_name']) != len(self._tool_inputs['FASTQ']) / 2:
            raise InvalidToolInputError(
                "Number of sample names does not correspond to the number of read files.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        reads_folder = self.__create_input_directory()
        self.__build_command(reads_folder)
        self._execute_command()
        self.__set_output()
        self.__analyze_metrics_file()

    def __get_reference(self) -> Path:
        """
        Returns the reference FASTA file, a symbolic link is created so the pipeline can index it.
        :return: Path to reference FASTA
        """
        link_path = self.folder / self._tool_inputs['FASTA'][0].path.name
        logger.info(f"Creating symlink for reference FASTA file: {link_path}")
        if not os.path.exists(link_path):
            os.symlink(self._tool_inputs['FASTA'][0].path, link_path)
        return link_path

    def __create_input_directory(self) -> Path:
        """
        Creates the input for the pipeline.
        :return: Reads input folder
        """
        dir_reads = Path(self._folder) / 'reads'
        if not dir_reads.exists():
            dir_reads.mkdir()
        for i, io_sample in zip(range(0, len(self._tool_inputs['FASTQ']), 2), self._tool_inputs['VAL_name']):
            forward_reads = Path(self._tool_inputs['FASTQ'][i].path)
            reverse_reads = Path(self._tool_inputs['FASTQ'][i + 1].path)
            sample_name_valid = fileutils.make_valid(io_sample.value)
            logger.info(f"Adding sample '{sample_name_valid}' as input")
            dir_sample = dir_reads / sample_name_valid
            if dir_sample.is_dir():
                shutil.rmtree(str(dir_sample))
            dir_sample.mkdir()
            if not fileutils.is_gzipped(forward_reads):
                (dir_sample / f"{sample_name_valid}_1.fastq").symlink_to(forward_reads)
                (dir_sample / f"{sample_name_valid}_2.fastq").symlink_to(reverse_reads)
            else:
                fileutils.gzip_extract(forward_reads, dir_sample / f"{sample_name_valid}_1.fastq")
                fileutils.gzip_extract(reverse_reads, dir_sample / f"{sample_name_valid}_2.fastq")
        return dir_reads

    def __build_command(self, input_dir: Path) -> None:
        """
        Builds the command.
        :param input_dir: Input directory
        :return: None
        """
        path_config = str(files('camel').joinpath('resources/tools/cfsan/snppipeline.conf'))
        command_parts = [
            'export CLASSPATH=$GATK_JAR:$PICARD_JAR:$CLASSPATH;',
            self._tool_command,
            f'-s {input_dir}',
            f'--conf {path_config}',
            str(self.__get_reference()),
            *self._build_options()
        ]
        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the output.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self._folder / 'snpma.fasta')]
        self._tool_outputs['FASTA_Preserved'] = [ToolIOFile(self._folder / 'snpma_preserved.fasta')]

        for key in ['VCF_Cons', 'VCF_Cons_preserved', 'VCF', 'VCF_preserved', 'BAM', 'Sample_names']:
            self._tool_outputs[key] = []
        with open(os.path.join(self._folder, 'sampleDirectories.txt')) as handle:
            for line in handle.readlines():
                sample_directory = Path(line.strip())
                self._tool_outputs['VCF_Cons'].append(ToolIOFile(sample_directory / 'consensus.vcf'))
                self._tool_outputs['VCF_Cons_preserved'].append(
                    ToolIOFile(sample_directory / 'consensus_preserved.vcf'))
                self._tool_outputs['VCF'].append(ToolIOFile(sample_directory / 'var.flt.vcf'))
                self._tool_outputs['VCF_preserved'].append(
                    ToolIOFile(sample_directory / 'var.flt_preserved.vcf'))
                self._tool_outputs['BAM'].append(ToolIOFile(sample_directory / 'reads.sorted.bam'))
                self._tool_outputs['Sample_names'].append(ToolIOValue(sample_directory.name))

    def __analyze_metrics_file(self) -> None:
        """
        Analyzes the metrics file and adds the information to the informs.
        :return: None
        """
        metrics_file = os.path.join(self._folder, 'metrics.tsv')
        if not os.path.isfile(metrics_file):
            raise FileNotFoundError("No metrics file generated.")
        with open(metrics_file) as handle:
            lines = handle.readlines()
        column_names = lines[0].split('\t')
        for line in lines[1:]:
            line_parts = line.split('\t')
            sample_name = line_parts[0].replace('"', '')
            stats = {column_names[i]: line_parts[i] for i in range(0, len(line_parts))}
            self._informs[sample_name] = stats

    def _check_command_output(self, command: Command) -> None:
        """
        Checks the command output.
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
