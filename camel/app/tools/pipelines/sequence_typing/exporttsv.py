from pathlib import Path

import pandas as pd
from camelcore.app.io.tooliofile import ToolIOFile

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.sequencetyping.typinghitbase import TypingHitBase


class ExportTSV(Tool):
    """
    Exports the typing output in TSV format.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Export TSV', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'HITS' not in self._tool_inputs:
            raise InvalidToolInputError('HITS input is required')
        super()._check_input()

    def _export_tsv_file(self, path_out: Path, include_hashed: bool) -> None:
        """
        Exports a TSV file.
        :param path_out: Output path
        :param include_hashed: If True, hashed alleles are included
        :return: None
        """
        hits: list[TypingHitBase] = [h.value for h in self._tool_inputs['HITS']]
        data_hits = pd.DataFrame(
            data=[h.to_table_row() for h in hits],
            columns=hits[0].table_column_names()
        )
        data_hits.to_csv(path_out, sep='\t', index=False)

    def _execute_tool(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        if len(self._tool_inputs['HITS']) == 0:
            logger.info('No hits found, not exporting TSV file.')
            return

        # Export the TSV file
        path_out = self._folder / self.get_param_value('output_filename')
        self._export_tsv_file(path_out, include_hashed=False)
        self._tool_outputs['TSV'] = [ToolIOFile(path_out)]

        # Export the TSV file with hashes
        if not any([h.value.is_new_allele() for h in self._tool_inputs['HITS']]):
            logger.debug('No novel alleles found, not exporting hashed TSV file.')
            return
        path_out_hash = path_out.parent / f'{path_out.stem}-hashes.tsv'
        self._export_tsv_file(path_out_hash, include_hashed=True)
