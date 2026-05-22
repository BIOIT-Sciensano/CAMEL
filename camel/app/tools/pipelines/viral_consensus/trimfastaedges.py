import logging

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.tool import Tool


class TrimFastaEdges(Tool):
    """
    Trims the edges of the input FASTA sequences when they are not covered.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('Trim FASTA edges', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError('FASTA input is required')
        if 'BED' not in self._tool_inputs:
            raise InvalidToolInputError('BED input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input BED file
        data_regions = pd.read_table(
            self._tool_inputs['BED'][0].path, names=['seq_id', 'start', 'stop'], keep_default_na=False, na_values='-')
        data_regions['size'] = data_regions['stop'] - data_regions['start']
        logging.info(f'{len(data_regions)} region(s) parsed')

        # Parse genome data
        with open(self._tool_inputs['FASTA'][0].path) as handle:
            seq_record_by_id = {s.id: s for s in SeqIO.parse(handle, 'fasta')}
        logging.info(f'{len(seq_record_by_id)} sequences parsed')
        data_genome = self.__extract_genome_data(seq_record_by_id)
        len_by_seq_id = {seq_id: length for seq_id, length in zip(data_genome['seq_id'], data_genome['length'])}
        self._informs['len_by_seq_id_in'] = len_by_seq_id

        # Mark regions at the edges
        if len(data_regions) > 0:
            data_regions['is_edge'] = data_regions.apply(
                lambda x: TrimFastaEdges._is_edge(x, len_by_seq_id[x['seq_id']]), axis=1)

        # Track stats
        self._informs['nb_clipped'] = 0
        self._informs['nb_masked'] = 0

        # Trim edges and replace other gaps by N's
        seq_records_out = []
        for id_, seq_record in seq_record_by_id.items():
            sequence = str(seq_record.seq)
            for region in data_regions.to_dict('records'):
                if region['seq_id'] != id_:
                    continue
                # Replace nucleotides in low-depth regions by Ns
                seq_before = sequence[:region['start']]
                seq_after = sequence[region['stop']:]
                sequence = seq_before + ('N' * region['size']) + seq_after

                # Update stats
                key_informs = 'nb_clipped' if region['is_edge'] else 'nb_masked'
                self._informs[key_informs] += region['size']

            # Strip Ns at the start and end
            seq_records_out.append(SeqRecord(Seq(sequence.strip('N')), id=id_, description=''))

        # Create output FASTA file
        path_fasta_out = self.folder / 'reference_masked.fasta'
        with path_fasta_out.open('w') as handle:
            SeqIO.write(seq_records_out, handle, 'fasta')
        self._informs['len_by_seq_id_out'] = {s.id: len(s) for s in seq_records_out}
        self._tool_outputs['FASTA'] = [ToolIOFile(path_fasta_out)]

    @staticmethod
    def _is_edge(record: pd.Series, seq_len: int) -> bool:
        """
        Returns true if the input region is at the edge of the input sequence (i.e., start or end).
        :param record: Input record
        :param seq_len: Sequence length
        :return: True if at edge, False otherwise
        """
        if record['start'] == 0:
            return True
        elif record['stop'] == seq_len:
            return True
        return False

    def __extract_genome_data(self, seq_record_by_id: dict[str, SeqIO.SeqRecord]) -> pd.DataFrame:
        """
        Creates a genome file needed by bedtools.
        :return: Path to genome file
        """
        records_out = []
        for seq_id, seq_record in sorted(seq_record_by_id.items()):
            records_out.append({'seq_id': seq_id, 'length': len(seq_record)})
        return pd.DataFrame(records_out)
