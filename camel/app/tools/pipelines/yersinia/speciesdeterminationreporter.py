import pandas as pd
from pathlib import Path

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool
from camel.app.error import InvalidToolInputError


class SpeciesDeterminationReporter(Tool):
    """
    This class is used to generate an output report for the species determination.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Yersinia: species determination reporter', '1.0')
        self._section = HtmlReportSection('Species determination')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if self._input_informs['analysis']['found_match']:
            self._section.add_header("Best match", level=3)
            table_df = pd.read_csv(self._tool_inputs['TSV_analysis'][0].path, sep="\t", header=0).fillna('NA')
            species, lineage, biotype, serotype, match, threshold = self._input_informs['analysis']['best_match'].values()
            self._section.add_table(data=[['Species:', f'<i>{species}</i>'],['Lineage:', f'{lineage}'],
                                     ['Biotype:', f'{biotype}'],['Serotype:', f'{serotype}']], column_names=None,
                                    table_attributes=[('class', 'information')])
            self._section.add_horizontal_line()
            self._section.add_header("All matches", level=3)
            self.__add_table_detected_species(table_df)
        else:
            self._section.add_paragraph('No <i>Yersinia</i> species found.')
        self._section.add_paragraph("Species and lineage designations as defined by <a href=\"https://doi.org/10.1099%2Fmgen.0.000301\">Savin et al</a>.")
        relative_path = Path('species_determination') / self._tool_inputs['TSV_analysis'][0].path.name
        self._section.add_link_to_file('Download results (TSV)', relative_path)
        self._section.add_file(self._tool_inputs['TSV_analysis'][0].path, relative_path)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def _check_input(self) -> None:
        """
        Checks if the provided tool input is valid.
        :return: None
        """
        if 'TSV_analysis' not in self._tool_inputs:
            raise InvalidToolInputError("TSV file with analysis results is required.")
        super()._check_input()

    def __add_table_detected_species(self, table_df) -> None:
        """
        Adds the table with the detected species.
        :return: None
        """
        table_html = []
        for _, r in table_df.iterrows():
            match = r['match']
            threshold = r['threshold']
            row = [f"<i>{r['species']}</i>", r['lineage'], r['biotype'], r['serotype'], f'{match:.2f}', f'{threshold:.2f}']
            table_html.append(row)
        header = ['Species', 'Lineage', 'Biotype', 'Serotype', 'Proportion of loci matched', 'Threshold for species and lineage']
        self._section.add_table(table_html, header, [('class', 'data')])
