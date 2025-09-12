from pathlib import Path

import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error import InvalidToolInputError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class NextcladeReporter(Tool):
    """
    Creates an HTML report for the Nextclade analysis.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Nextclade reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'CSV' not in self._tool_inputs:
            raise InvalidToolInputError('Nextclade CSV input is required')
        if 'DB' not in self._tool_inputs:
            raise InvalidToolInputError('Nextclade DB input is required')
        if 'nextclade' not in self._input_informs:
            raise InvalidToolInputError('Nextclade informs are required')
        if len(self._input_informs['nextclade']['results']) != 1:
            logger.warning(f'{self.name} only supports a single nextclade result')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('Nextclade', subtitle=self._input_informs['nextclade']['_name'])
        informs_in = self._input_informs['nextclade']['results'][0]

        # Quality control
        section.add_header('Quality control', 4)
        section.add_table([
            ['Overall QC score:', f"{int(informs_in['qc.overallScore']):,}"],
            ['Overall QC status:', informs_in['qc.overallStatus']],
            ['Mixed sites:', informs_in['qc.mixedSites.status']],
            ['Private mutations:', informs_in['qc.privateMutations.status']],
            ['Frameshifts:', informs_in['qc.frameShifts.status']],
            ['Stop codons:', informs_in['qc.stopCodons.status']],
        ], None, [('class', 'information')])

        # Clade
        # Retrieve column from the DB folder
        column_mapping = pd.read_table(self._tool_inputs['DB'][0].path / 'report_clade_cols.tsv')
        section.add_header('Clade', 4)
        section.add_table([
            [f'{name}:', informs_in[key]] for key, name in zip(column_mapping['key'], column_mapping['name'])],
            None, [('class', 'information')])

        # Mutations
        section.add_header('Mutations', 4)
        section.add_table([
            ['Nucleotide'] + [informs_in[k] for k in (
                'totalSubstitutions', 'totalInsertions', 'totalDeletions')],
            ['Amino-acid'] + [informs_in[k] for k in (
                'totalAminoacidSubstitutions', 'totalAminoacidInsertions', 'totalAminoacidDeletions')],
        ], ['Level', 'Substitutions', 'Insertions', 'Deletions'], [('class', 'data')])
        section.add_paragraph('A complete overview of all detected mutations is available in the CSV file below.')

        # Overview link
        section.add_header('Overview', 4)
        section.add_paragraph('<b>Note:</b> Nextclade only analyzes the HA segment.')
        relative_path = Path('nextclade', f"nextclade_{self._parameters['name'].value}.tsv")
        section.add_link_to_file('Download (CSV)', relative_path)
        section.add_file(self._tool_inputs['CSV'][0].path, relative_path)

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
