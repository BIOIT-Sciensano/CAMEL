from typing import Optional

from camel.app.toolkits.pubmlst.pubmlstparser import PubMLSTParser, PubMLSTParsingError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class ResistanceMetadataExtractor(Tool):
    """
    This tool is used to extract metadata for resistance genes in the Neisseria pipeline.
    """

    def __init__(self) -> None:
        """
        Initializes the metadata extractor.
                :return: None
        """
        super().__init__('Neisseria: resistance metadata extractor', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        loci = self._parameters['loci'].value.split(', ')
        table_data = []
        for h in self._tool_inputs['hits']:
            if h.value.locus not in loci:
                continue
            allele_url = self.__get_allele_url(h.value.locus, h.value.allele_id)
            if allele_url is None:
                logger.warning("Cannot determine valid url for: {}_{}".format(h.value.locus, h.value.allele_id))
                continue
            try:
                key, value = PubMLSTParser.parse_linked_data(allele_url)
            except PubMLSTParsingError as err:
                logger.warning(f"Cannot parse: {allele_url}: {err}")
                continue
            table_data.append(['{} ({}_{}):'.format(key, h.value.locus, h.value.allele_id), value])

        section = self._tool_inputs['VAL_HTML'][0].value
        section.add_table(table_data, table_attributes=[('class', 'information')])
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(section)]

    def __get_allele_url(self, locus: str, allele_id: str) -> Optional[str]:
        """
        Retrieves the allele url for the given allele.
        :param locus: Locus
        :param allele_id: Allele id
        :return: Allele url
        """
        locus_informs = self._input_informs['scheme']['loci'].metadata_by_locus_name[locus.lower()]
        if allele_id == '-' or allele_id is None:
            return None
        return locus_informs['allele_page_url'].format(allele_id=allele_id)
