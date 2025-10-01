from pathlib import Path

import pydotplus

from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class RdCsbReporter(Tool):
    """
    This class reports the output for the CSB / RD1 based species identification for Mycobacterium.
    """

    TITLE = '<i>csb</i> & regions of difference (RD1, RD9)'

    VISUALIZATION_PARAMS = {
        'color_present': '#ccebc5',
        'color_absent': '#fbb4ae',
        'shape_width': 1.4,
        'shape_height': 0.5
    }

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: RD / csb reporter', '0.1')
        self._sub_directory = Path('rd_csb')

    def _check_input(self) -> None:
        """
        Checks if the tool input is valid.
        :return: None
        """
        if 'HITS' not in self._tool_inputs:
            raise InvalidToolInputError("Gene detection hits input (HITS) is required")
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("FASTA input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection(RdCsbReporter.TITLE, 3)
        self.__set_informs(self._tool_inputs['HITS'])
        self.__add_visualization()
        self.__add_output_tables(self._tool_inputs['HITS'])
        self.__add_warning()
        self.__add_database(self._tool_inputs['FASTA'][0].path)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def __set_informs(self, hits: list) -> None:
        """
        Sets the informs for this tool.
        :param hits: Gene detection hits
        """
        self._informs['loci_detected'] = [h.value.locus for h in hits]
        if 'csb' in self._informs['loci_detected']:
            species = 'M. tuberculosis' if 'RD9' in self._informs['loci_detected'] else 'M. africanum'
        else:
            species = 'Other M. bovis' if 'RD1' in self._informs['loci_detected'] else 'M. bovis BCG'
        self._informs['species'] = species

    @staticmethod
    def __create_node(text: str, italic: bool = False, is_present: bool = False) -> pydotplus.Node:
        """
        Creates a node for the visualization.
        :param text: Node text
        :param italic: If True, text is shown in italic
        :param is_present: Determines the color for the node
        :return: Node
        """
        return pydotplus.Node(
            text,
            width=RdCsbReporter.VISUALIZATION_PARAMS['shape_width'],
            height=RdCsbReporter.VISUALIZATION_PARAMS['shape_height'],
            fontname='times' if not italic else 'times italic',
            style='filled',
            fillcolor=RdCsbReporter.VISUALIZATION_PARAMS['color_present' if is_present else 'color_absent'])

    def __add_visualization(self) -> None:
        """
        Creates a visualization of the decision tree.
        :return: None
        """
        graph = pydotplus.Dot(graph_type='graph')

        # First level
        node_csb = pydotplus.Node('csb', fontname='times italic', shape='diamond')
        graph.add_node(node_csb)

        # Second level
        node_bovine = RdCsbReporter.__create_node('Bovine', False, 'csb' not in self._informs['loci_detected'])
        graph.add_node(node_bovine)
        graph.add_edge(pydotplus.Edge(node_csb, node_bovine, label='not present', fontname='times'))

        node_non_bovine = RdCsbReporter.__create_node('Non-bovine', False, 'csb' in self._informs['loci_detected'])
        graph.add_node(node_non_bovine)
        graph.add_edge(pydotplus.Edge(node_csb, node_non_bovine, label='present', fontname='times'))

        # Second level (non-bovine)
        node_rd9 = pydotplus.Node('RD9', shape='diamond', fontname='times')
        graph.add_node(node_rd9)
        graph.add_edge(pydotplus.Edge(node_non_bovine, node_rd9))

        node_afr = RdCsbReporter.__create_node('M. africanum', True, self._informs['species'] == 'M. africanum')
        graph.add_node(node_afr)
        graph.add_edge(pydotplus.Edge(node_rd9, node_afr, label='not present', fontname='times'))

        node_tbc = RdCsbReporter.__create_node('M. tuberculosis', True, self._informs['species'] == 'M. tuberculosis')
        graph.add_node(node_tbc)
        graph.add_edge(pydotplus.Edge(node_rd9, node_tbc, label='present', fontname='times'))

        # Second level (bovine)
        node_rd1 = pydotplus.Node('RD1', shape='diamond', fontname='times')
        graph.add_node(node_rd1)

        edge_rd1 = pydotplus.Edge(node_bovine, node_rd1)
        graph.add_edge(edge_rd1)

        node_bcg = RdCsbReporter.__create_node('M. bovis BCG', True, self._informs['species'] == 'M. bovis BCG')
        graph.add_node(node_bcg)
        graph.add_edge(pydotplus.Edge(node_rd1, node_bcg, label='not present', fontname='times'))

        node_bovis = RdCsbReporter.__create_node('M. bovis', True, self._informs['species'] == 'Other M. bovis')
        graph.add_node(node_bovis)
        graph.add_edge(pydotplus.Edge(node_rd1, node_bovis, label='present', fontname='times'))

        output_path = self._folder / 'graph-csb-rd.png'
        graph.write(str(output_path), format='png')

        # Add to the report
        self._report_section.add_header('Interpretation', 4)
        relative_path = self._sub_directory / output_path.name
        self._report_section.add_file(output_path, relative_path)
        img = HtmlElement('img', attributes=[(
            'class', 'bordered'), ('src', str(relative_path)), ('alt', 'visualization')])
        self._report_section.add_html_object(img)

    def __add_output_tables(self, hits: list) -> None:
        """
        Adds the output table which shows the presence/absence of RD1 and csb.
        :param hits: Gene detection hits
        :return: None
        """
        hits_by_locus = {h.value.locus: h.value for h in hits}
        self._report_section.add_header("Detected sequences", 4)
        header = ['Locus', 'Present']
        for name in self._input_informs['columns'][1:]:
            header.append(name)
        table_data = []
        for locus in ['csb', 'RD1', 'RD9']:
            locus_text = '<i>csb</i>' if locus == 'csb' else locus
            if locus in hits_by_locus:
                hit = hits_by_locus[locus]
                table_data.append([HtmlTableCell(locus_text, hit.color), HtmlTableCell('yes', hit.color)] +
                                  hits_by_locus[locus].to_html_row(self._report_section, self._sub_directory, )[1:])
            else:
                table_data.append([HtmlTableCell(locus_text, 'red'), HtmlTableCell('no', 'red')] +
                                  [HtmlTableCell('-', 'red')] * (len(header) - 2))

        # Remove accession column
        index = header.index('Accession')
        for row in table_data:
            row.pop(index)
        header.pop(index)

        self._report_section.add_table(table_data, header, [('class', 'data')])
        self._report_section.add_alert(
            f"Detection for this DB is always done using 'KMA' regardless of pipeline setting.", 'info')

    def __add_database(self, fasta_path: Path) -> None:
        """
        Adds a download link for the database to the report.
        :param fasta_path: FASTA path
        :return: None
        """
        relative_path = self._sub_directory / fasta_path.name
        self._report_section.add_file(fasta_path, relative_path)
        self._report_section.add_link_to_file('Database (FASTA)', relative_path)

    def __add_warning(self) -> None:
        """
        Adds a warning to the report when pseudo-reads were used as input.
        :return: None
        """
        if 'pseudo_reads' in self._parameters:
            self._report_section.add_warning_message("The tool is executed on simulated reads.")

    @staticmethod
    def generate_empty_section() -> HtmlReportSection:
        """
        Returns an empty report for when this analysis is disabled.
        :return: Report section
        """
        section = HtmlReportSection(RdCsbReporter.TITLE)
        section.add_paragraph('Analysis disabled')
        return section
