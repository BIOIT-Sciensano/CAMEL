from pathlib import Path

from camel.app.core.utils import fileutils
from camel.app.toolkits.export.tsvexporter import TsvExporter
from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlexpandablediv import HtmlExpandableDiv
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.toolkits.mycobacterium import assay51snputils
from camel.app.toolkits.mycobacterium.assay51snputils import SNPPosition
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class Assay51SnpReporter(Tool):
    """
    This class reports the detected species based on a set of 51 SNPs.
    The SNPs and the methodology are described in:
    https://jcm.asm.org/content/52/6/1962.full
    """

    TITLE = '51 SNP-based'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: 51SNP reporter', '0.1')
        self._section = HtmlReportSection(Assay51SnpReporter.TITLE)
        self._sub_folder = Path('51snp')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__add_overview_table()
        positions_by_name = self._input_informs['detection']['snp_positions_by_name']
        self.__add_gyrb_table(positions_by_name)
        self.__add_genetic_group_table(positions_by_name)
        self.__add_scg_table(positions_by_name)
        self.__add_snps_table(positions_by_name)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_overview_table(self) -> None:
        """
        Adds a table with an overview of the assay.
        :return: None
        """
        table_data = [
            ['MTBC positive control:', self._input_informs['detection']['mtbc_pos_control']],
            ['<i>gyrB</i> species differentiation:', '{} ({})'.format(
                self._input_informs['detection']['gyrB_group'], self._input_informs['detection']['gyrB_species'])],
            ['Genetic group:', self._input_informs['detection']['genetic_group']],
            ['Best matching SCG:', '{} (ST{}) - {}/45 SNPs ({:.2f}% match)'.format(
                self._input_informs['detection']['scg_profile'].scg,
                self._input_informs['detection']['scg_profile'].st,
                self._input_informs['detection']['scg_nb_snps_matched'],
                100 * self._input_informs['detection']['scg_nb_snps_matched'] / 45)]
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_gyrb_table(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Adds the table with the result of the gyrB assay.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        table_data = []
        for profile in assay51snputils.GYRB_PROFILES:
            color = 'yellow' if profile['group'] == self._input_informs['detection']['gyrB_group'] else None
            row = [HtmlTableCell(profile['group'], color=color)]
            for key in ['SNP02', 'SNP03', 'SNP04']:
                row.append(HtmlTableCell(profile[key], color=snp_positions_by_name[key].get_color(profile[key])))
            row.append(HtmlTableCell(profile['species'], color=color))
            table_data.append(row)
        header = ['Group', '<i>gyrB</i> 1450 (SNP02)', '<i>gyrB</i> 756 (SNP03)', '<i>gyrB</i> 675 (SNP04)', 'Species']
        self._section.add_header('<i>gyrB</i> species differentiation', 3)
        self._section.add_table(table_data, header, [('class', 'data')])

    def __add_genetic_group_table(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Adds the table with the result of the genetic group assay.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        table_data = []
        for genetic_group in assay51snputils.GENETIC_GROUPS:
            color = 'yellow' if genetic_group['name'] == self._input_informs['detection']['genetic_group'] else None
            row = [HtmlTableCell(genetic_group['name'], color)]
            for key in ['SNP05', 'SNP06']:
                color = snp_positions_by_name[key].get_color(genetic_group[key])
                row.append(HtmlTableCell(genetic_group[key], color=color))
            table_data.append(row)
        header = ['Genetic group', '<i>katG</i> 463 (SNP05)', '<i>gyrA</i> 95 (SNP06)']
        self._section.add_header('Genetic group', 3)
        self._section.add_table(table_data, header, [('class', 'data')])

    def __add_scg_table(self, snp_positions_by_name: dict[str, SNPPosition], nb_columns: int = 6) -> None:
        """
        Adds a table with the SNPs supporting the detected SCG.
        :param snp_positions_by_name: SNP positions by name
        :param nb_columns: Nb. of columns
        :return: None
        """
        table_data = [[HtmlElement('th', 'Matching SNPs', [('colspan', nb_columns)])]]
        row = []
        for i, snp in enumerate(self._input_informs['detection']['scg_profile'].snps):
            position = snp_positions_by_name[f'SNP{i+7:02d}']
            row.append(HtmlTableCell(position.name, color=position.get_color(snp)))
            if len(row) > 5:
                table_data.append(row)
                row = []
        table_data.append(row)
        self._section.add_header('Best matching SNP cluster group (SCG)', 3)
        self._section.add_table(table_data, table_attributes=[('class', 'data')])

    def __add_snps_table(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Adds a table with an overview of all SNPs.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        table_data = []
        for name, pos in sorted(snp_positions_by_name.items()):
            row = [name, pos.pos, pos.ref]
            if pos.vcf_record is not None:
                row.append('{}{}'.format(str(pos.vcf_record.ALT[0]), '*' if pos.vcf_filt_record is None else ''))
            else:
                row.append('-')
            table_data.append(row)
        header = ['Name', 'Position', 'Ref. (H37Rv)', 'Alt.']
        div = HtmlExpandableDiv('51snps', 'SNPs')
        div.add_table(table_data, header, [('class', 'data')])
        self._section.add_header('All 51 SNPs', 3)
        self._section.add_html_object(div)

        table_path = self._folder / 'all_snps-{}.tsv'.format(
            fileutils.make_valid(self._tool_inputs['VAL_Sample'][0].value))
        TsvExporter.export(table_data, header, table_path)
        relative_path = self._sub_folder / table_path.name
        self._section.add_file(table_path, relative_path)
        self._section.add_link_to_file('Download (TSV)', relative_path)
        self._section.add_paragraph(
            "This table contains the SNPs compared to the reference genome, SNPs indicated with a '*' did not pass "
            "filtering.")
