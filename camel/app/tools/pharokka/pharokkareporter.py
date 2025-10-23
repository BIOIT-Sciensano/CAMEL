from pathlib import Path

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.reports.htmltableformatter import HtmlTableFormatter, FormatEntry
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class PharokkaReporter(Tool):
    """
    Parses Pharokka output files and generates an HTML report section for the final report.
    """

    TITLE = 'Pharokka'
    SUB_DIR = 'pharokka'

    COLS_HITS: list[FormatEntry] = [
        {'key': 'contig', 'title': 'Contig'},
        {'key': 'gene', 'title': 'Gene'},
        {'key': 'hit', 'title': 'Hit'},
        {'key': 'alnScore', 'title': 'Alignment score'},
        {'key': 'seqIdentity', 'title': 'Sequence identity (%)', 'fmt': lambda x: f'{x * 100:.2f}'},
        {'key': 'start', 'title': 'Start'},
        {'key': 'stop', 'title': 'Stop'},
        {'key': 'frame', 'title': 'Frame'},
    ]

    COLS_GENERAL: list[FormatEntry] = [
        {'key': 'contig', 'title': 'Contig'},
        {'key': 'Accession', 'title': 'Accession'},
        {'key': 'Description', 'title': 'Description'},
        {'key': 'Classification', 'title': 'Classification'},
        {'key': 'Modification_Date', 'title': 'Modif. date'},
        {'key': 'mash_distance', 'title': 'mash distance', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'mash_pval', 'title': 'mash p-value', 'fmt': lambda x: f'{x:.4f}'},
        {'key': 'mash_matching_hashes', 'title': 'mash matching hashes'},
    ]

    COLS_GENOMIC_PROPERTIES: list[FormatEntry] = [
        {'key': 'contig', 'title': 'Contig'},
        {'key': 'Genome_Length_(bp)', 'title': 'Length', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'molGC_(%)', 'title': '%GC', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'Coding_Capacity_(%)', 'title': 'Coding capacity (%)', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'Low_Coding_Capacity_Warning', 'title': 'Low coding capacity warning',
         'fmt': lambda x: 'Yes' if not pd.isna(x) else 'No'},
        {'key': 'Positive_Strand_(%)', 'title': 'Pos. strand (%)', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'Negative_Strand_(%)', 'title': 'Neg. strand (%)', 'fmt': HtmlTableFormatter.FLOAT_FMT},
        {'key': 'Number_CDS', 'title': 'Nb. of CDS', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'tRNAs', 'title': 'Nb. of tRNAs', 'fmt': HtmlTableFormatter.INT_FMT},
        {'key': 'Jumbophage', 'title': 'Low coding capacity warning', 'fmt': lambda x: 'Yes' if x is True else 'No'},
    ]

    COLS_HOST: list[FormatEntry] = [
        {'key': 'contig', 'title': 'Contig'},
        {'key': 'Host', 'title': 'Host'},
        {'key': 'Isolation_Host_(beware_inconsistent_and_nonsense_values)', 'title': 'Isolation host*'}
    ]

    COLS_TAXONOMY: list[FormatEntry] = [
        {'key': 'contig', 'title': 'Contig'},
        {'key': 'Lowest_Taxa', 'title': 'Lowest taxa'},
        {'key': 'Realm', 'title': 'Realm'},
        {'key': 'Kingdom', 'title': 'Kingdom'},
        {'key': 'Phylum', 'title': 'Phylum'},
        {'key': 'Class', 'title': 'Class'},
        {'key': 'Order', 'title': 'Order'},
        {'key': 'Family', 'title': 'Family'},
        {'key': 'Sub-family', 'title': 'Sub-family'},
        {'key': 'Genus', 'title': 'Genus'},
        {'key': 'Baltimore_Group', 'title': 'Baltimore group'},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Pharokka Reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'TSV_STATS' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka output is required ('pharokka_cds_functions.tsv')")
        if 'TSV_CARD' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka output is required ('CARD hits')")
        if 'TSV_VFDB' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka output is required ('VFDB hits')")
        if 'TSV_INPHARED' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka output is required ('pharokka_top_hits_mash_inphared.tsv')")
        if 'GBK' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka output is required ('Genbank file')")
        if 'PNG' not in self._tool_inputs:
            raise InvalidToolInputError("Pharokka Multiplotter output is required ('PNG')")
        if 'pharokka' not in self._input_informs:
            raise InvalidToolInputError("Pharokka informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection(PharokkaReporter.TITLE, subtitle=self._input_informs['pharokka']['_name'])

        # Summary metrics
        section.add_header('Annotation and stats', 3)
        self._add_download_table(section)
        section.add_horizontal_line()

        # Genomic maps
        section.add_header('Genomic map(s)', 3)
        self._add_download_plots_table(section)
        section.add_horizontal_line()

        # Antibiotic sensitivity
        if 'show_card' in self._parameters:
            section.add_header('Antibiotic susceptibility (CARD database)', 3)
            self.__add_antibiotic_sensitivity(section)
            section.add_horizontal_line()

        # Virulence factors
        if 'show_vfdb' in self._parameters:
            section.add_header('Virulence factor identification (VFDB database)', 3)
            self.__add_virulence_factors(section)
            section.add_horizontal_line()

        # Phage identification
        if 'show_inphared' in self._parameters:
            section.add_header('Phage identification (INPHARED database)', 3)
            self.__add_identification(section)
            section.add_horizontal_line()

        # Store the output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def _add_download_table(self, section: HtmlReportSection) -> None:
        """
        Adds links to the summary metrics and the genbank file.
        :param section: Report output section
        :return: None
        """
        # Add files
        rel_path_stats = Path(PharokkaReporter.SUB_DIR, 'annotation_stats.tsv')
        section.add_file(self._tool_inputs['TSV_STATS'][0].path, rel_path_stats)
        rel_path_gbk = Path(PharokkaReporter.SUB_DIR, self._tool_inputs['GBK'][0].path.name)
        section.add_file(self._tool_inputs['GBK'][0].path, rel_path_gbk)

        # Add table
        download_table = [
            ['Genbank file', HtmlTableCell('Download (GBK)', link=str(rel_path_gbk))],
            ['Annotation stats', HtmlTableCell('Download (TSV)', link=str(rel_path_stats))]
            ]
        section.add_table(download_table, ['File', 'Download'], [('class', 'data')])

    def _add_download_plots_table(self, section: HtmlReportSection) -> None:
        """
        Adds links to the genomic map(s).
        :param section: Report output section
        :return: None
        """
        # Add plot(s)
        rel_paths_png = {}
        for file in self._tool_inputs['PNG']:
            rel_paths_png[file] = Path(PharokkaReporter.SUB_DIR, '/pharokka_plots/', file.path.name)
        for file, rel_path in rel_paths_png.items():
            section.add_file(file.path, rel_path)

        gen_map_table = []
        for i in range(len(rel_paths_png)):
            rel_path_png = list(rel_paths_png.values())[i]
            if i == 0:
                new_map = ['Genomic map', HtmlTableCell('Download (PNG)', link=str(rel_path_png))]
            else:
                new_name = 'Genomic map ' + str(i+1)
                new_map = [new_name, HtmlTableCell('Download (PNG)', link=str(rel_path_png))]
            gen_map_table.append(new_map)

        # Add table
        section.add_table(gen_map_table, ['File', 'Download'], [('class', 'data')])

    def __add_antibiotic_sensitivity(self, section: HtmlReportSection) -> None:
        """
        Adds the table with CARD AMR hits.
        :param section: Report output section
        :return: None
        """
        data_hits = pd.read_table(self._tool_inputs['TSV_CARD'][0].path)
        if data_hits.empty:
            section.add_paragraph('No AMR genes were detected.')
        else:
            # Update column names
            data_hits.columns = [col.split('_')[-1] for col in data_hits.columns]
            header = [c['title'] for c in PharokkaReporter.COLS_HITS]
            section.add_table(HtmlTableFormatter.format_table_data(
                data_hits, PharokkaReporter.COLS_HITS), header, [('class', 'data')])

    def __add_virulence_factors(self, section: HtmlReportSection) -> None:
        """
        Adds the table with VFDB hits.
        :param section: Report output section
        :return: None
        """
        data_hits = pd.read_table(self._tool_inputs['TSV_VFDB'][0].path)
        logger.info(list(data_hits.columns))
        if data_hits.empty:
            section.add_paragraph('No virulence genes were detected.')
        else:
            header = [c['title'] for c in PharokkaReporter.COLS_HITS]
            section.add_table(HtmlTableFormatter.format_table_data(
                data_hits, PharokkaReporter.COLS_HITS), header, [('class', 'data')])

    def __add_identification(self, section: HtmlReportSection) -> None:
        """
        Adds the table with the INPHARED information.
        :param section: Report output section
        :return: None
        """
        data = pd.read_table(self._tool_inputs['TSV_INPHARED'][0].path)

        # Add output table (split by category)
        section.add_header('General information', 4)
        section.add_table(HtmlTableFormatter.format_table_data(
            data, PharokkaReporter.COLS_GENERAL), [c['title'] for c in PharokkaReporter.COLS_GENERAL],
            [('class', 'data')])

        section.add_header('Genomic properties', 4)
        section.add_table(HtmlTableFormatter.format_table_data(
            data, PharokkaReporter.COLS_GENOMIC_PROPERTIES),
            [c['title'] for c in PharokkaReporter.COLS_GENOMIC_PROPERTIES], [('class', 'data')])

        section.add_header('Host', 4)
        section.add_table(HtmlTableFormatter.format_table_data(
            data, PharokkaReporter.COLS_HOST), [c['title'] for c in PharokkaReporter.COLS_HOST], [('class', 'data')])
        section.add_paragraph('(*) Beware of inconsistent and nonsense values.')

        section.add_header('Taxonomy', 4)
        section.add_table(HtmlTableFormatter.format_table_data(
            data, PharokkaReporter.COLS_TAXONOMY), [c['title'] for c in PharokkaReporter.COLS_TAXONOMY],
            [('class', 'data')])

        section.add_header('Download', 4)
        path_rel = Path('pharokka', 'inphared.tsv')
        section.add_link_to_file('Download (TSV)', path_rel)
        section.add_file(self._tool_inputs['TSV_INPHARED'][0].path, path_rel)
