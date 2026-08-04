from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.reports.htmltablecell import HtmlTableCell

from camel.app.core.tool import Tool


class SerogroupDeterminationReporter(Tool):
    """
    This class is used to generate an output report for the serogroup determination.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Neisseria: serogroup determination reporter', '0.1')
        self._section = HtmlReportSection('Serogroup determination (Legacy)')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        detected_serogroup = self._input_informs['analysis']['detected_serogroup']
        self._section.add_paragraph(f'Detected serogroup: <b>{detected_serogroup}</b>')
        self.__add_table_detected_loci()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_table_detected_loci(self) -> None:
        """
        Adds the table with the detected loci.
        :return: None
        """
        table_data = []
        for s in self._input_informs['analysis']['serogroups_sorted']:
            row = [s['name'], '{:.2f}'.format(100 * s['fraction_detected'])]
            for i in range(0, 7):
                try:
                    locus_name, color = s['color_per_hit'][i]
                    row.append(HtmlTableCell(locus_name, color))
                except IndexError:
                    row.append('-')
            table_data.append(row)
        header = ['Serogroup', 'Genes detected (%)'] + [''] * 7
        self._section.add_table(table_data, header, [('class', 'data')])
