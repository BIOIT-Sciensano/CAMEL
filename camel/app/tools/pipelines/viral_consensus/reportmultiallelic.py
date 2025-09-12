from pathlib import Path

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ReporterMultiAllelic(Tool):
    """
    Creates an output report for the multi-allelic site calling.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Reporter: Multi-allelic site calling', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Multi-allelic sites input file is required (TSV)')
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('Updated consensus sequencing input is required (FASTA)')
        if 'calling' not in self._input_informs:
            raise InvalidToolInputError('Calling informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        section = HtmlReportSection('Multi-allelic sites')
        try:
            data_multi_allelic_sites = pd.read_table(self._tool_inputs['TSV'][0].path)
        except pd.errors.EmptyDataError:
            data_multi_allelic_sites = []

        # Add results section
        section.add_header('Results', 4)
        section.add_table([
            ['Nb. multi-allelic sites', f'{len(data_multi_allelic_sites):,}'],
        ], ['Metric', 'Value'], [('class', 'data')])

        # Add output files section
        section.add_header('Output files', 4)
        relative_path_tsv = Path('multi_allelic', f"overview_multi_allelic_sites-{self._parameters['name'].value}.tsv")
        section.add_file(self._tool_inputs['TSV'][0].path, relative_path_tsv)
        relative_path_fasta = Path('multi_allelic', f"{self._parameters['name'].value}-ambiguous.fasta")
        section.add_file(self._tool_inputs['FASTA'][0].path, relative_path_fasta)

        section.add_table([
            ['Overview', HtmlTableCell('Download (TSV)', link=str(relative_path_tsv))],
            ['Updated consensus sequence', HtmlTableCell('Download (FASTA)', link=str(relative_path_fasta))]],
            ['File', 'Download'], [('class', 'data')]
        )

        # Add the information section
        section.add_header('Additional information', 4)
        section.add_paragraph(
            'Multi-allelic sites are identified by creating a pileup using <i>bcftools</i>. Positions where the '
            'minority allele has an allele frequency of at least <b>{:.2f}</b> are classified as multi-allelic and '
            'encoded in the updated consensus sequence according to the IUPAC codes. Minimum depth: <b>{}x</b>.'.format(
                self._input_informs['calling']['min_freq_minor_allele'], self._input_informs['calling']['min_dp'])
        )

        # Store output
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]
