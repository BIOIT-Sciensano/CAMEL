from pathlib import Path
from typing import Any

from camel.app.core.reports.htmlreport import HtmlReport
from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe.fastqinput import FastqInput
from camel.app.scriptutils.basescript.scriptinput import ScriptInput

from camel.app.scriptutils.inputhelper.inputhelperbase import InputHelperBase, TAssemblyOpts, TTrimmingOpts


class FastaHelper(InputHelperBase):
    """
    Helper for FASTA input.
    """

    @staticmethod
    def opts_from_cli(opts: dict[str, Any]) -> tuple[TAssemblyOpts | None, TTrimmingOpts | None, int]:
        """
        Parses the options from the command line.
        :return: Assembly options, trimming options, number of threads
        """
        return None, None, 1

    def prepare_fasta_input(self, script_in: ScriptInput, report: HtmlReport) -> Path:
        """
        Prepares the FASTA input.
        :param script_in: Script input
        :param report: HTML report
        :return: Path to FASTA file
        """
        if not script_in.type_ == model.InputType.FASTA:
            raise ValueError(f"Invalid input type: {script_in.type_}")
        dir_symlinks = self._dir / 'symlinks'
        dir_symlinks.mkdir(exist_ok=True, parents=True)
        path_symlink = dir_symlinks / f'{self._name}.fasta'
        path_symlink.symlink_to(script_in.fasta)
        return path_symlink

    def prepare_fastq_input(self, fq_in: FastqInput, report: HtmlReport) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param fq_in: FASTQ input
        :param report: HTML report
        :return: FASTQ input
        """
        raise NotImplementedError("FASTA input does not support FASTQ input")
