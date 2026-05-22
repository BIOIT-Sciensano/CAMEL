from pathlib import Path

from camelcore.app.io.tooliovalue import ToolIOValue
from camelcore.app.reports.htmlreportsection import HtmlReportSection
from camelcore.app.utils import fileutils

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class HtmlReporterAssembly(Tool):
    """
    Tool to create HTML reports for the Assembly.
    """

    def __init__(self) -> None:
        """
        Initialize this tool.
        :return: None
        """
        super().__init__('HTML Reporter', '0.1')
        self.__subfolder = Path('assembly')
        self._report_section = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection('Assembly', subtitle=self._input_informs['assembler'])
        self.__add_assembly_info()
        self.__add_assembly_download_link()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section, False)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTA_raw' not in self._tool_inputs:
            raise InvalidToolInputError("No assembly input found ('FASTA_raw')")
        if 'FASTA_filt' not in self._tool_inputs:
            raise InvalidToolInputError("No filtered assembly input found ('FASTA_filt')")
        if 'SAMPLE_NAME' not in self._tool_inputs:
            raise InvalidToolInputError("No sample name input found ('SAMPLE_NAME')")
        if 'quast' not in self._input_informs:
            raise InvalidToolInputError("Quast informs are required ('quast')")
        if 'assembler' not in self._input_informs:
            raise InvalidToolInputError("Assembler informs are required ('assembler')")
        super()._check_input()

    def __add_assembly_info(self) -> None:
        """
        Adds the assembly info.
        :return: None
        """
        quast_informs = self._input_informs['quast']
        table_data = [
            ('Assembler:', self._input_informs['assembler']['_name']),
            ('N50:', '{:,}'.format(int(quast_informs['contig']['N50']))),
            ('Number of contigs:', '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)']))),
            ('Total length:', '{:,}'.format(int(quast_informs['genome']['Total length'])))
        ]
        self._report_section.add_table(table_data, table_attributes=[('class', 'information')])

    def __add_assembly_download_link(self) -> None:
        """
        Adds a download link for the assembly.
        :return: None
        """
        sample_name_valid = fileutils.make_valid(self._tool_inputs['SAMPLE_NAME'][0].value)

        # Add filtered assembly
        relative_path = self.__subfolder / f'{sample_name_valid}_contigs.fasta'
        self._report_section.add_file(self._tool_inputs['FASTA_filt'][0].path, relative_path)
        self._report_section.add_link_to_file('Assembly (FASTA)', relative_path)

        # Add unfiltered assembly
        relative_path_raw = self.__subfolder / f'{sample_name_valid}_contigs_unfilt.fasta'
        self._report_section.add_file(self._tool_inputs['FASTA_raw'][0].path, relative_path_raw)
