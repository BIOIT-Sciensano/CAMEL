import pandas as pd
import yaml
from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class NextcladeSubTypeReporter(Tool):
    """
    Creates an HTML report for the Nextclade subtype determination.
    """

    COLUMNS = [
        {'key': 'identity', 'title': 'Identity', 'fmt': lambda x: f'{x * 100:.2f}%'},
        {'key': 'hashes', 'title': 'Matching hashes'},
        {'key': 'm_mult', 'title': 'Median-multiplicity (Cov.)'},
        {'key': 'q_id', 'title': 'Query', 'fmt': lambda x: x.replace('.fasta', '')},
    ]

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Nextclade subtype reporter', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidToolInputError('Nextclade TSV input is required')
        if 'mash' not in self._input_informs:
            raise InvalidToolInputError('mash informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('Subtype determination', subtitle=self._input_informs['mash']['_name'])
        if len(self._tool_inputs['TSV']) == 0:
            section.add_alert('Subtype determination failed', 'error')
            subtype = None
        else:
            data_in = pd.read_table(self._tool_inputs['TSV'][0].path, names=self._input_informs['mash']['cols'])
            if len(data_in) == 0:
                section.add_alert('Subtype determination failed', 'error')
                subtype = None
                self.informs['nextclade_dbs'] = {}
            else:
                # Extract the best matching subtype
                data_in.sort_values(by='identity', inplace=True, ascending=False)
                query = data_in.iloc[0]['q_id']
                subtype = query.split('-')[0]
                path_yml = self._tool_inputs['DB'][0].path / 'nextclade_dbs.yml'
                if not path_yml.exists():
                    raise FileNotFoundError(f'Nextclade database mapping not found: {path_yml}')
                with path_yml.open() as handle:
                    data = yaml.safe_load(handle)
                self.informs['nextclade_dbs'] = data.get(subtype, {})

                # Report output
                section.add_paragraph(f"Detected subtype: <b>{subtype}</b>")
                section.add_table(
                    [[row[col['key']] if col.get('fmt') is None else col['fmt'](row[col['key']]) for
                      col in NextcladeSubTypeReporter.COLUMNS] for row in data_in.to_dict('records')],
                    [x['title'] for x in NextcladeSubTypeReporter.COLUMNS], [('class', 'data')])
                section.add_alert(
                    'When the correct subtype is not automatically detected, the correct <i>nextclade</i> database can '
                    'be specified in the interface.', 'info')

        # Tool output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
        self._informs['subtype'] = subtype
