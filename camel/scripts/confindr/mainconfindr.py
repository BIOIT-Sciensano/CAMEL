#!/usr/bin/env python
import dataclasses
import json
from pathlib import Path
from typing import Any

import click
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.tools.confindr.confindr import ConFindr
from camel.app.tools.confindr.confindrreporter import ConFindrReporter
from camel.resources import DIR_CITATIONS


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Custom options for the ConFindr script.
    """

    db: Path = dataclasses.field(metadata={'help': 'Path to ConFindr database'})
    working_dir: Path = dataclasses.field(
        default=Path.cwd(), metadata={'help': 'Working directory'}
    )
    rmlst: bool = dataclasses.field(
        default=False,
        metadata={
            'help': 'Prefer using rMLST databases over core-gene derived databases'
        },
    )
    quality_cutoff: int = dataclasses.field(
        default=20, metadata={'help': 'Base quality cutoff'}
    )
    base_cutoff: int = dataclasses.field(
        default=3, metadata={'help': 'Number of bases  cutoff'}
    )
    base_percentage_cutoff: int = dataclasses.field(
        default=5, metadata={'help': 'Base percentage cutoff'}
    )
    min_matching_hashes: int = dataclasses.field(
        default=150, metadata={'help': 'Minimum number of matching KMA hashes'}
    )
    threads: int = 1


class MainConFindr(BaseScript[ScriptInput, ScriptOutput, Options]):
    """
    Main script for the ConFindr tool.
    """

    def __init__(
        self, script_in: ScriptInput, script_out: ScriptOutput, script_opts: Options
    ) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param script_out: Script output
        ;param script_opts: Script options
        :return: None
        """
        super().__init__(
            name='ConFindr',
            version='1.0',
            script_in=script_in,
            script_out=script_out,
            script_opts=script_opts,
        )

    def _execute(self) -> None:
        """
        Executes this script.
        :return: None
        """
        # Initialize the report
        report = reportutils.init_report(
            path_out=self._script_out.html,
            key=self._name,
            title='ConFindr',
            dir_out=self._script_out.dir,
        )
        report.add_html_object(
            reportutils.create_overview_section(
                version=self.version,
                dataset_name=self._script_in.name,
                input_file_str=self._script_in.input_str,
                input_type=self._script_in.type_.value,
            )
        )
        report.save()

        # Check if the database exists
        if not self._script_opts.db.exists():
            raise FileNotFoundError(f'DB not found: {self._script_opts.db}')

        # Run ConFindr
        input_dict = self.__prepare_input()
        confindr = ConFindr()
        confindr.update_parameters(
            databases=str(self._script_opts.db),
            quality_cutoff=self._script_opts.quality_cutoff,
            base_cutoff=self._script_opts.base_cutoff,
            base_fraction_cutoff=self._script_opts.base_percentage_cutoff / 100,
            min_matching_hashes=self._script_opts.min_matching_hashes,
            data_type={
                model.InputType.ONT: 'Nanopore',
                model.InputType.ILLUMINA: 'Illumina',
            }[self._script_in.type_],
            rmlst=self._script_opts.rmlst,
            threads=self._script_opts.threads,
        )
        if self._script_opts.rmlst:
            confindr.update_parameters(rmlst=True)
        confindr.add_input_files(input_dict)
        confindr.run(self._script_opts.working_dir.absolute())

        # Save informs to file (if specified)
        if self._script_out.json is not None:
            with self._script_out.json.open('w') as handle:
                json.dump(confindr.informs, handle, indent=2)
                logger.info(f'ConFindr informs saved to {self._script_out.json}')

        # Create output report
        confindr_reporter = ConFindrReporter()
        confindr_reporter.add_input_informs({'confindr': confindr.informs})
        confindr_reporter.update_parameters(input_type=self._script_in.type_.value)
        confindr_reporter.run(self._script_opts.working_dir.absolute())
        report.add_html_object(confindr_reporter.tool_outputs['HTML'][0].value)

        # Add citation and command
        report.add_html_object(
            reportutils.create_commands_section(
                [confindr.informs], self._script_opts.working_dir
            )
        )
        report.add_html_object(
            reportutils.create_citations_section(
                dir_=DIR_CITATIONS,
                keys_other=['Low_2019-confindr', 'Jolley_2012-rmlst'],
            )
        )
        report.save()
        logger.info(f'Report saved to: {self._script_out.html}')

    def __prepare_input(self) -> dict[str, Any]:
        """
        Prepares the input for the confindr tool.
        :return: Input dictionary
        """
        if self._script_in.fastq_pe is not None:
            return {'FASTQ_PE': [ToolIOFile(fq) for fq in self._script_in.fastq_pe]}
        else:
            return {'FASTQ_SE': [ToolIOFile(self._script_in.fastq_se)]}


@click.command(name='confindr', short_help='Wrapper for ConFindr')
@basescriptutils.add_input_opts(
    supported=[model.InputType.ILLUMINA, model.InputType.ONT]
)
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Entry point for the common interface.
    :param kwargs: Command line arguments
    :return: None
    """
    script = MainConFindr(
        script_in=basescriptutils.parse_script_input(kwargs),
        script_out=basescriptutils.parse_script_output(kwargs),
        script_opts=Options(**cliutils.from_kwargs(Options, kwargs)),
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
