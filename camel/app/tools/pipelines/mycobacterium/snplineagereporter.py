from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class SNPLineageReporter(Tool):
    """
    This tool is used to generate output report sections.
    """

    TITLE = 'SNP lineage'

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: SNP lineage reporter', '0.1')
        self._section = HtmlReportSection('SNP lineage')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'detection' not in self._input_informs:
            raise InvalidToolInputError("Linaege detection informs are required ('detection')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__add_lineage_overview_table()
        self.__add_variants_table()
        self.__add_db_info()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_lineage_overview_table(self) -> None:
        """
        Adds a table with an overview of the lineages.
        :return: None
        """
        self._section.add_header('Predicted lineage(s)', 4)
        header = ['Predicted lineage', 'Full name', 'Main spoligotype', 'RDS', '# supporting SNPs']
        table_data = []
        for level, data in self._input_informs['detection']['detected_lineage_by_level'].items():
            if data is None:
                continue
            lineage = data['lineage']
            table_data.append([lineage.id_, lineage.name, lineage.main_spoligo, lineage.rd_type, data['count']])
        self._section.add_table(table_data, header, [('class', 'data')])

    def __add_variants_table(self) -> None:
        """
        Adds the table with the lineage associated SNPs.
        :return: None
        """
        self._section.add_header('Lineage supporting SNPs', 4)
        data = []
        header = ['Position', 'Ref', 'Sample', 'Lineage', 'Passes filtering']
        for snp in self._input_informs['detection']['detected_snps']:
            if snp.lineage.id_ in ('lineage4', 'lineage4.9'):
                data.append([snp.start, snp.ref, snp.ref, snp.lineage.id_, '-'])
            else:
                data.append([snp.start, snp.ref, snp.alt, snp.lineage.id_, 'Yes' if snp.passes_filtering else 'No'])
        self._section.add_table(data, header, [('class', 'data')])
        self._section.add_paragraph("For lineage 4 and sublineage 4.9 the nucleotide at the SNP positions should not "
                                    "contain a SNP, because the reference genome belongs to lineage 4.9.")

    def __add_db_info(self) -> None:
        """
        Adds the database information.
        """
        self._section.add_header('Database information', level=4)
        self._section.add_table([
            ['Origin', 'TBDB (TBProfiler)'],
            ['Last updated', self._input_informs['detection']['db_version']]
        ], ['Field', 'Value'], [('class', 'data')])
