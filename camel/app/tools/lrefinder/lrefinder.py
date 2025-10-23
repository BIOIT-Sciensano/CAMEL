import re
from typing import Any

import bs4

from camel.app.core.command import Command
from camel.app.core.utils import toolutils
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class LREFinder(Tool):
    """
    Tool for detection of the 23S rRNA mutations, and the optrA, cfr, cfr(B) and poxtA genes, encoding linezolid
    resistance in Enterococci from whole genome sequences.
    Reference: https://cge.cbs.dtu.dk/services/LRE-Finder/
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('LRE-Finder', '20200812')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in ('FASTQ_SE', 'FASTQ_PE')):
            raise InvalidToolInputError('FASTQ_(SE/PE) input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        output_path = self.folder / 'lrefinder_out.tsv'
        if 'FASTQ_SE' in self._tool_inputs:
            fq_in = f"-i {self._tool_inputs['FASTQ_SE'][0].path}"
        else:
            fq_in = f"-ipe {' '.join(str(io.path) for io in self._tool_inputs['FASTQ_PE'])}"
        # create command itself
        self._command.command = ' '.join([
                self._tool_command,
                fq_in,
                '-t_db', '$LREFINDER_DB',
                '-o', str(output_path)
        ])
        self._execute_command()
        self._informs.update(self.__parse_html_output(self._command.stdout))

    # noinspection PyTypeChecker
    def __parse_html_output(self, html_code: str) -> dict[str, Any]:
        """
        Parses the HTML output reported by the tool.
        :param html_code: HTML code
        :return: Parsed information
        """
        info = {}
        soup = bs4.BeautifulSoup(html_code, 'html.parser')
        info['species'] = soup.find('b', string='Species identified:').findNext('table').find('th').text
        info['genes'] = LREFinder.__parse_html_table(soup.find('b', string='Genes Identified:').findNext('table'))
        info['mutations'] = LREFinder.__parse_html_table(
            soup.find('b', string='Identified mutations in 23s:').findNext('table'))
        return info

    @staticmethod
    def __parse_html_table(html_table: bs4.ResultSet) -> list[dict[str, Any]]:
        """
        Parses an HTML table.
        :param html_table: HTML table
        :return: Parsed information
        """
        # noinspection PyUnresolvedReferences
        rows = html_table.findAll('tr')
        header = [th.text.replace('_', ' ') for th in rows[0].findAll('th')]
        header = [re.sub(r'\[(.*)]', r'(\g<1>)', value) for value in header]
        table_data = [{
            col: td.text for col, td in zip(header, row.findAll('td'))} for row in rows[1:] if
            len(row.findAll('td')) > 0]

        # Convert numeric values to float
        for info_dict in table_data:
            for key in info_dict.keys():
                try:
                    info_dict[key] = float(info_dict[key])
                except ValueError:
                    continue
        return table_data

    def _check_command_output(self, command: Command) -> None:
        """
        Checks if the command was executed successfully.
        :param command: Command to check
        :return: None
        """
        toolutils.check_tool_execution(self, command, exit_code=0)
