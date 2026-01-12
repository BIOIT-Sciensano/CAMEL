import pandas as pd

from camel.app.core.reports.htmlelement import HtmlElement
from camel.app.core.reports.htmlexpandablediv import HtmlExpandableDiv
from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.mobsuite.mobreconreporter import MOBReconReporter
from camel.app.core.tool import Tool


class GenomicContext(Tool):
    """
    Tool to link detected genes to their genomic context (e.g. plasmid / chromosome).
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Genomic context', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'dbs' not in self._input_informs:
            raise InvalidToolInputError('Database informs input is required')
        super()._check_input()

    def _get_plasmid_status(self, contig_name: str, plasmids: list[str]) -> list[HtmlTableCell]:
        """
        Returns the plasmid status for the given contig.
        :param contig_name: Name of the contig
        :param plasmids: List of plasmids
        :return: List of cells for the output table
        """
        contig_status = self._input_informs['mob_recon']['contig_report'].get(contig_name)
        cells = [HtmlTableCell('X', color='green') if contig_status is None else HtmlTableCell('-')]
        for plasmid in plasmids:
            cells.append(HtmlTableCell('X', color='green') if plasmid == contig_status else HtmlTableCell('-'))
        return cells

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse mob_recon input data
        plasmids = sorted(list(x for x in set(
            self._input_informs['mob_recon']['contig_report'].values()) if x is not None))

        # Create the output report
        section = HtmlReportSection('Genomic context', level=3)

        # Check if the detection method is BLAST
        if self._parameters['detection_method'].value != 'blast':
            logger.warning('Genomic context can only be predicted with blast as detection method')
            section.add_paragraph('Predicting genomic context is only performed when the detection method is blast.')
        elif len(self._input_informs['dbs']) == 0:
            section.add_paragraph('No compatible databases selected.')
        else:
            if self._parameters['input_type'] == 'illumina':
                section.add_warning_message(
                    'Predicting genomic context based solely on short-read data is error-prone and should only be '
                    'considered as an indication.')

            # Create rows
            for db in self._input_informs['dbs']:
                section.add_header(db['title'], 3)

                # No hits found
                if db['idx'] is None:
                    section.add_paragraph('No hits found.')
                    continue

                # Parse hits
                data_hits = pd.read_table(self._tool_inputs['TSV'][db['idx']].path)
                if len(data_hits) == 0:
                    section.add_paragraph('No hits found.')
                    continue

                # Add table with hits
                if len(data_hits) > 10:
                    div = HtmlExpandableDiv(f"genomic_context-{db['key']}", f'{len(data_hits)} rows.')
                else:
                    div = HtmlElement('div')

                header = ['Key', 'Chromosome', *[MOBReconReporter.format_plasmid_id(p) for p in plasmids]]
                div.add_table([
                    [f"<i>{row[db['gene']]}</i>",
                     *self._get_plasmid_status(row[db['contig']], plasmids)] for row in data_hits.to_dict('records')
                ], header, [('class', 'data')])
                section.add_html_object(div)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
