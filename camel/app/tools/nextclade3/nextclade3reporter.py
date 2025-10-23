from pathlib import Path

import pandas as pd

from camel.app.core.reports.htmlreportsection import HtmlReportSection
from camel.app.core.reports.htmltablecell import HtmlTableCell
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.core.tool import Tool


class Nextclade3Reporter(Tool):
    """
    Reporter for Nextclade 3.
    """

    DIR = 'nextclade'

    KEYS_QC = [
        ('qc.overallStatus', 'Overall QC status'),
        ('qc.mixedSites.status', 'Mixed sites'),
        ('qc.privateMutations.status', 'Private mutations'),
        ('qc.frameShifts.status', 'Frameshifts'),
        ('qc.stopCodons.status', 'Stop codons')
    ]

    KEYS_MUTS = [
        ('totalSubstitutions', '# substitutions (Nucl.)'),
        ('totalInsertions', '# insertions (Nucl.)'),
        ('totalDeletions', '# deletions (Nucl.)'),
        ('totalAminoacidSubstitutions', '# substitutions (AA)'),
        ('totalAminoacidInsertions', '# insertions (AA)'),
        ('totalAminoacidDeletions', '# deletions (AA)')
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Nextclade3 reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Nextclade TSV input is required')
        if 'nextclade' not in self._input_informs:
            raise InvalidToolInputError('Nextclade informs are required')
        super()._check_input()

    @staticmethod
    def format_qc_cell(value: str) -> HtmlTableCell:
        """
        Formats the value into a QC cell.
        :param value: Input value
        :return: Formatted table cell
        """
        if pd.isna(value):
            return HtmlTableCell('-', color='grey')
        elif value == 'good':
            return HtmlTableCell('good', color='green')
        elif value == 'mediocre':
            return HtmlTableCell('mediocre', color='yellow')
        elif value == 'bad':
            return HtmlTableCell('bad', color='red')
        raise ValueError(f'Invalid QC value: {value}')

    def __get_output_name(self, segment: str) -> str:
        """
        Returns the output name for the Nextclade TSV output.
        :param segment: Segment
        :return: Output name
        """
        if 'name' in self._parameters:
            name = self._parameters['name'].value
            return f'nextclade-{name}_{segment}.tsv'
        return f'nextclade-{segment}.tsv'

    def __parse_tsv_input(self) -> pd.DataFrame:
        """
        Parses the TSV input files.
        :return: Parsed data
        """
        records_out = []
        for path_tsv in [x.path for x in self._tool_inputs['TSV']]:
            data = pd.read_table(path_tsv)
            data['segment'] = path_tsv.parent.name
            records_out.extend(data.to_dict('records'))
        return pd.DataFrame(records_out)

    def __add_table_qc(self, section: HtmlReportSection, data_nextclade: pd.DataFrame) -> None:
        """
        Adds the table with the QC overview.
        :param section: HTML report section
        :param data_nextclade: Nextclade data
        :return: None
        """
        section.add_header('Quality control', 4)
        header = ['Segment', 'Overall QC score', *[title for _, title in Nextclade3Reporter.KEYS_QC]]
        section.add_table([[
            row['segment'].upper() if 'capitalize_segment_names' in self._parameters else row['segment'],
            f"{row['qc.overallScore']:.2f}",
            *[Nextclade3Reporter.format_qc_cell(row[k]) for k, _ in Nextclade3Reporter.KEYS_QC]
        ] for row in data_nextclade.to_dict('records')], header, [('class', 'data')])

    def __add_table_mutations(self, section: HtmlReportSection, data_nextclade: pd.DataFrame) -> None:
        """
        Adds the table with the mutation overview.
        :param section: HTML report section
        :param data_nextclade: Nextclade data
        :return: None
        """
        section.add_header('Mutations', 4)
        header = ['Segment', *[title for _, title in Nextclade3Reporter.KEYS_MUTS]]
        section.add_table([[
            row['segment'].upper() if 'capitalize_segment_names' in self._parameters else row['segment'],
            *[row[k] for k, _ in Nextclade3Reporter.KEYS_MUTS]
        ] for row in data_nextclade.fillna('n/a').to_dict('records')], header, [('class', 'data')])
        section.add_paragraph('A complete overview of all detected mutations is available in the TSV output.')

    def __add_table_metadata(self, section: HtmlReportSection, data_nextclade: pd.DataFrame) -> None:
        """
        Adds the metadata table.
        :param section: HTML report section
        :param data_nextclade: Nextclade data
        :return: None
        """
        section.add_header('Metadata', 4)

        # Determine which columns need to be added
        columns = {}
        for informs in self._input_informs['nextclade']:
            cols_metadata = informs['db'].get('metadata_columns')
            if cols_metadata is None:
                continue
            for row in cols_metadata:
                if row['key'] in columns:
                    continue
                columns[row['key']] = row['name']

        # Add table
        header = ['Segment', *[title for _, title in columns.items()]]
        section.add_table([[
            row['segment'].upper() if 'capitalize_segment_names' in self._parameters else row['segment'],
            *[row[k] for k, _ in columns.items()]
        ] for row in data_nextclade.fillna('n/a').to_dict('records')], header, [('class', 'data')])


    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if len(self._input_informs['nextclade']) > 0:
            section = HtmlReportSection('Nextclade', subtitle=self._input_informs['nextclade'][0]['_name'])
            data_nextclade = self.__parse_tsv_input()
            self.__add_table_qc(section, data_nextclade)
            self.__add_table_mutations(section, data_nextclade)
            self.__add_table_metadata(section, data_nextclade)

            # Add download links
            section.add_header('Downloads', 4)
            header = ['Segment', 'Download (TSV)']
            table_data = []
            for i, segment in enumerate(data_nextclade['segment'].unique()):
                relative_path = Path(self.DIR, self.__get_output_name(segment))
                section.add_file(self._tool_inputs['TSV'][i].path, relative_path)
                segment_name = segment.upper() if 'capitalize_segment_names' in self._parameters else segment
                table_data.append([segment_name, HtmlTableCell('Download (TSV)', link=str(relative_path))])
            section.add_table(table_data, header, [('class', 'data')])

            # Add database information
            section.add_header('Database information', 4)
            header = ['Segment', 'Version', 'Reference']
            section.add_table([[
                segment.upper() if 'capitalize_segment_names' in self._parameters else segment,
                informs['db']['version'],
                informs['db']['reference']
            ] for segment, informs in zip(data_nextclade['segment'], self._input_informs['nextclade'])],
                header, [('class', 'data')])
        else:
            section = HtmlReportSection('Nextclade')
            section.add_paragraph('No matching nextclade databases found.')

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
