#!/usr/bin/env python
import dataclasses
import shutil
from pathlib import Path

import click
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.utils import fileutils, reportutils

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging
from camel.app.scriptutils.basescript import basescriptutils
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.basescript.fastainput import FastaInput
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput
from camel.app.scriptutils.model import BaseOptions
from camel.app.toolkits.genedetection.dbhelper import DBHelper
from camel.app.toolkits.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.tools.cdhit.cdhitest import Cluster


@dataclasses.dataclass(frozen=True)
class Options(BaseOptions):
    """
    Options for the make DB script.
    """

    working_dir: Path = dataclasses.field(
        default=Path.cwd(), metadata={'help': 'Working directory'}
    )
    threads: int = dataclasses.field(
        default=4, metadata={'help': 'Number of threads to use', 'show_default': True}
    )
    identity_cutoff: int = dataclasses.field(
        default=80, metadata={'help': 'Identity cutoff for clustering'}
    )


class MainMakeGeneDetectionDB(BaseScript[FastaInput, ScriptOutput, Options]):
    """
    Creates databases for the gene detection tool.
    """

    def __init__(self, in_: FastaInput, out: ScriptOutput, opts: Options) -> None:
        """
        Initializes the main script.
        :param in_: Script input
        :param out: Script output
        :param opts: Options
        :return: None
        """
        super().__init__(
            name='Gene detection: Make DB',
            version='1.0.0',
            script_in=in_,
            script_out=out,
            script_opts=opts,
        )
        self._db_name = fileutils.make_valid(self._script_in.name)
        self._helper = DBHelper(self._db_name, self._script_opts.working_dir)
        self._clusters: list[Cluster] | None = None
        self._new_name_by_header = None

    def _execute(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        if not self._script_out.dir.exists():
            self._script_out.dir.mkdir(parents=True)
        input_fasta = self._helper.standardize_fasta_headers(self._script_in.fasta)
        self.__export_blast_db(input_fasta, self._script_out.dir)
        self._clusters = self.__cluster_fasta(input_fasta)
        self._helper.export_metadata(self._db_name, self._script_out.dir)
        self.__export_report()

    def __cluster_fasta(self, input_fasta: Path) -> list[Cluster]:
        """
        Clusters the input FASTA file.
        :param input_fasta: Input FASTA file
        :return: List of clusters
        """
        dir_clustering = self._helper.get_working_subdir('clustering')
        fasta_seq_headers = dir_clustering / 'seq_headers.fasta'
        self._new_name_by_header = self._helper.convert_fasta_headers_to_seq(
            input_fasta, fasta_seq_headers
        )
        return self._helper.get_clusters_form_fasta(
            fasta_seq_headers,
            self._script_opts.identity_cutoff,
            self._script_opts.threads,
        )

    def __export_blast_db(self, input_fasta: Path, output_dir: Path) -> None:
        """
        Creates and exports a gene detection BLAST database from the given FASTA file.
        :param input_fasta: Input FASTA file
        :param output_dir: Output directory
        :return: None
        """
        # Create file
        dir_indexing = self._helper.get_working_subdir('index_blast')
        new_path = dir_indexing / input_fasta.name
        shutil.copyfile(str(input_fasta), str(new_path))

        # Index
        self._helper.index_samtools_faidx(new_path, dir_indexing)
        self._helper.index_blast(new_path, dir_indexing)

        # Export files
        for f in dir_indexing.iterdir():
            shutil.copyfile(str(f), str(output_dir / f.name))

    def __export_report(self) -> None:
        """
        Creates a report with some info on the database.
        :return: None
        """
        self._report = reportutils.init_report(
            self._script_out.html,
            'gene_detection',
            'Gene detection database',
            self._script_out.dir,
        )
        self._report.add_html_object(self.__create_db_info_section())
        self._report.add_html_object(self.__create_clusters_section())
        self._report.add_html_object(
            reportutils.create_commands_section(
                self._helper.informs, self._script_opts.working_dir
            )
        )
        self._report.save()

    def __create_db_info_section(self) -> HtmlReportSection:
        """
        Creates the report section with the database info.
        :return: HTML report section
        """
        section_db_info = HtmlReportSection('Database info')
        section_db_info.add_table(
            [
                ['Name:', self._db_name],
                ['Size:', sum(len(c.seq_ids) for c in self._clusters)],
                ['Nb. clusters:', len(self._clusters)],
                ['Clustering cutoff: ', f'{self._script_opts.identity_cutoff}%'],
            ],
            table_attributes=[('class', 'information')],
        )
        return section_db_info

    def __create_clusters_section(self) -> HtmlReportSection:
        """
        Creates the report section with the cluster information.
        :return: HTML report section
        """
        section_clusters = HtmlReportSection('Clusters')
        allele_by_seq_id = {
            s: GeneDetectionUtils.parse_header(h)[1]['allele']
            for s, h in self._new_name_by_header.items()
        }
        table_data = [
            [
                cluster.name,
                len(cluster.seq_ids),
                ', '.join([allele_by_seq_id[s] for s in cluster.seq_ids]),
            ]
            for cluster in self._clusters
        ]
        section_clusters.add_table(
            table_data, ['Cluster', 'Size', 'Sequence ids'], [('class', 'data')]
        )
        return section_clusters


@click.command(
    name='gene_detection_create_db',
    short_help='Creates DBs for the gene detection script',
)
@cliutils.add_click_options_from_dataclass(FastaInput)
@basescriptutils.add_output_opts
@cliutils.add_click_options_from_dataclass(Options)
def main(**kwargs) -> None:
    """
    Wrapper for NCBI AMRFinder+.
    """
    script = MainMakeGeneDetectionDB(
        in_=FastaInput(**cliutils.from_kwargs(FastaInput, kwargs)),
        out=basescriptutils.parse_script_output(kwargs),
        opts=Options(**cliutils.from_kwargs(Options, kwargs)),
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
