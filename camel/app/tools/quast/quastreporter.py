from pathlib import Path
from typing import Any

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool
from camel.app.core.utils import fileutils


class QuastReporter(Tool):
    """
    Creates an output report for QUAST with additional modules.
    """

    def __init__(self) -> None:
        """
        Initialize this tool.
        :return: None
        """
        super().__init__('QUAST reporter', '5.2.0')

    def _check_input(self) -> None:
        """
        Checks whether required quast TSV input is available
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required.")
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError("TSV input is required.")
        if 'HTML' not in self._tool_inputs:
            raise InvalidToolInputError("HTML input is required.")
        if 'DIR' not in self._tool_inputs:
            raise InvalidToolInputError("DIR input is required.")
        if 'quast' not in self._input_informs:
            raise InvalidToolInputError("QUAST informs are required.")
        if 'assembler' not in self._input_informs:
            raise InvalidToolInputError("Assembler informs are required.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Execute this tool.
        :return: None
        """
        subtitle = ', '.join([self._input_informs['quast']['_name_full'], self._input_informs['busco']['_name_full']])
        section = HtmlReportSection('Assembly quality', subtitle=subtitle)

        # Parse QUAST input file
        with self._tool_inputs['TSV'][0].path.open() as handle:
            data_quast = {parts[0]: parts[1] for parts in [line.strip().split('\t') for line in handle.readlines()]}

        # Basic statistics
        self.__add_section_basic_stats(section, data_quast)
        section.add_horizontal_line()

        # Completeness
        self.__add_section_completeness(section, data_quast)
        section.add_horizontal_line()

        # Coverage
        self.__add_section_coverage(section, data_quast)
        section.add_horizontal_line()

        # Output files
        self.__add_section_downloads(section)
        self.__copy_quast_dir(section)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]

    def __add_section_basic_stats(self, section: HtmlReportSection, data_quast: dict[str, Any]) -> None:
        """
        Adds the section with the basic statistics.
        :param section: Report section
        :param data_quast: QUAST data
        :return: None
        """
        section.add_header('Basic statistics', 3)
        section.add_table([
            ['Assembler:', self._input_informs['assembler']],
            ['Nb. of contigs:', f"{int(data_quast['# contigs']):,}"],
            ['Total length:', f"{int(data_quast['Total length']):,}"],
            ['N50:', f"{int(data_quast['N50']):,}"]
        ], None, [('class', 'information')])
        section.add_paragraph(
            '<b>Note:</b> Contigs shorter than 500 bp are not considered for the calculation of these metrics.')

    def __add_section_completeness(self, section: HtmlReportSection, data_quast: dict[str, Any]) -> None:
        """
        Adds the section with the completeness statistics.
        :param section: Report section
        :param data_quast: QUAST data
        :return: None
        """
        busco_stats = self._input_informs['busco']['results']['results'] if 'busco' in self._input_informs else None
        section.add_header('Completeness', 3)
        section.add_table([
            ['Reference length:', f"{int(data_quast['Reference length']):,}" if
                'Reference length' in data_quast else 'n/a'],
            ['Genome fraction:', f"{float(data_quast['Genome fraction (%)']):.2f}%" if
                'Genome fraction (%)' in data_quast else 'n/a'],
            ['Duplication ratio:', f"{float(data_quast['Duplication ratio']):.2f}" if
                'Duplication ratio' in data_quast else 'n/a'],
            ['Complete BUSCO:', f"{busco_stats.get('Complete'):.2f}%" if busco_stats else 'n/a'],
            ['Partial BUSCO:', f"{busco_stats.get('Fragmented'):.2f}%" if busco_stats else 'n/a'],
        ], None, [('class', 'information')])

    def __add_section_coverage(self, section: HtmlReportSection, data_quast: dict[str, Any]) -> None:
        """
        Adds the section with the coverage statistics.
        :param section: Report section
        :param data_quast: QUAST data
        :return: None
        """
        section.add_header('Coverage', 3)
        ref_genome = self._input_informs['quast']['ref'].replace('.fasta', '') if (
                'ref' in self._input_informs['quast']) else 'n/a'
        section.add_paragraph(f"Reference genome (RefSeq accession): {ref_genome}")
        section.add_table([
            [
                'Assembly',
                data_quast.get('Avg. coverage depth', 'n/a'),
                f"{float(data_quast['Coverage >= 1x (%)']):.2f}%" if 'Coverage >= 1x (%)' in data_quast else 'n/a'
            ], [
                'Reference',
                data_quast.get('Reference avg. coverage depth', 'n/a'),
                f"{float(data_quast['Reference coverage >= 1x (%)']):.2f}%" if
                'Reference coverage >= 1x (%)' in data_quast else 'n/a'
            ]
        ], ['Category', 'Avg. coverage', 'Positions covered >1x'], [('class', 'data')])

    def __add_section_downloads(self, section: HtmlReportSection) -> None:
        """
        Adds the section with the downloads.
        :param section: Report section
        :return: None
        """
        section.add_header('Downloads', 3)
        name_sanitized = fileutils.make_valid(self._parameters['name'].value)

        # QUAST report
        relative_path_html = Path('assembly', 'quast', f'quast_{name_sanitized}.html')
        section.add_file(self._tool_inputs['HTML'][0].path, relative_path_html)

        # TSV output
        relative_path_tsv = Path('assembly', 'quast', f'quast_{name_sanitized}.tsv')
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path_tsv)

        # Assembly
        relative_path_fasta = Path('assembly', f'{name_sanitized}.fasta')
        section.add_file(self._tool_inputs['FASTA'][0].path, relative_path_fasta)

        # Add table
        section.add_table([
            ['Assembly', HtmlTableCell('Download (FASTA)', link=str(relative_path_fasta))],
            ['QUAST report', HtmlTableCell('Download (HTML)', link=str(relative_path_html))],
        ], ['File', 'Download'], [('class', 'data')])

    def __copy_quast_dir(self, section: HtmlReportSection) -> None:
        """
        Copies the content of the QUAST directory.
        :param section: Output section
        :return: None
        """
        kept_extensions = ['pdf', 'txt', 'tsv', 'html']
        for p in self._tool_inputs['DIR'][0].path.rglob('*'):
            if not any(p.name.endswith(x) for x in kept_extensions):
                continue
            if p.name in ('report.html', 'report.tsv'):
                # Skip report & TSV file because they are already added
                continue
            path_rel = p.relative_to(self._tool_inputs['DIR'][0].path)
            section.add_file(p, Path('assembly', 'quast', path_rel))
