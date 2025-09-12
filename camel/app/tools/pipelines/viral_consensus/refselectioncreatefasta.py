import json
import logging
import math
from io import StringIO
from pathlib import Path
from typing import Union, Any

import pandas as pd
from Bio import SeqIO
from pysam.libcvcf import defaultdict

from camel.app.command.command import Command
from camel.app.components import toolutils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class RefSelection(Tool):
    """
    Parses the mash TSV output and selects the best matching reference sequence.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('viral consensus: reference selection', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        toolutils.check_input(self, keys_required=['TSV', 'DB'])
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input
        db_info = self._parse_database_info(self._tool_inputs['DB'][0].path)

        # Parse the length information
        len_by_segment = RefSelection._parse_seq_lengths(self._tool_inputs['DB'][0].path)

        # Parse the mash output
        mash_out_by_segment = RefSelection._parse_mash_output([p.path for p in self._tool_inputs['TSV']], len_by_segment)

        # Select the best reference for each segment
        ref_by_segment = {}
        for seg, mash_out in mash_out_by_segment.items():
            try:
                best_hit = mash_out.iloc[0]
                ref_by_segment[seg] = best_hit
            except IndexError:
                ref_by_segment[seg] = None

        # Eliminate mutually exclusive elements
        if 'mutually_exclusive_segments' in db_info:
            self._filter_mutually_exclusive_segments(ref_by_segment, db_info['mutually_exclusive_segments'])

        # Create merged FASTA & TSV files
        self._tool_outputs['FASTA'] = [ToolIOFile(self.__create_fasta(ref_by_segment, self._tool_inputs['DB'][0].path))]
        self._tool_outputs['TSV'] = [ToolIOFile(self._merge_mash_output(mash_out_by_segment))]

        # Save output in JSON format (selected references)
        path_out = self.folder / 'selected_refs.json'
        with path_out.open('w') as handle:
            json.dump({
                s: (r.to_dict() if r is not None else None) for s, r in ref_by_segment.items()}, handle, indent=2)
        self._tool_outputs['JSON'] = [ToolIOFile(path_out)]

    @staticmethod
    def calc_score(row: pd.Series) -> float:
        """
        Calculates the score for the mash output.
        It considers:
        - % identity
        - coverage (median multiplicity)
        - total hashes (otherwise shorter sequences are benefited)
        - percentage of hash matches
        :param row: Input row
        :return: Score
        """
        # noinspection PyTypeChecker
        return (
            (0.3 * row['identity']) +
            (0.2 * row['hashes_pct']) +
            (0.2 * math.log2(int(row['median_mult']) + 1)) +
            (0.2 * math.log2(int(row['length']) + 1))
         )

    @staticmethod
    def _parse_mash_output(mash_output: list[Path], len_by_segment: dict[str, dict[str, int]]) -> dict[str, pd.DataFrame]:
        """
        Parses the mash output for the input files.
        :param mash_output: mash output files
        :return: Parsed data by segment
        """
        mash_out_by_segment = {}
        for path_tsv in mash_output:
            segment = path_tsv.parents[1].name
            data_mash = pd.read_table(
                path_tsv, names=['identity', 'hashes', 'median_mult', 'p_val', 'ref_id', 'ref_comment'])
            mash_out_by_segment[segment] = data_mash
            data_mash['hashes_pct'] = data_mash['hashes'].apply(
                lambda x: 100 * int(x.split('/')[0]) / max(int(x.split('/')[1]), 1))
            data_mash['hashes_nb'] = data_mash['hashes'].apply(
                lambda x: int(x.split('/')[1]))
            data_mash['length'] = data_mash['ref_id'].map(len_by_segment[segment])
            data_mash['score_final'] = data_mash.apply(RefSelection.calc_score, axis=1)
            data_mash.sort_values(by=['score_final'], inplace=True, ascending=False)
            data_mash['ref_id_fmt'] = data_mash['ref_id'].apply(lambda x: x.split('-')[0].replace('_', ' '))
        return mash_out_by_segment

    @staticmethod
    def _parse_database_info(dir_db: Path) -> dict:
        """
        Parses the database information.
        :param dir_db: Database directory
        """
        path_metadata = dir_db / 'genome_info.json'
        with path_metadata.open() as handle:
            return json.load(handle)

    @staticmethod
    def _parse_seq_lengths(dir_db: Path) -> dict[str, dict[str, int]]:
        """
        Parse the lenghts of the sequences in the database.
        :param dir_db: Database directory
        :return: Sequence length by segment
        """
        len_by_segment = defaultdict(dict)
        for path_fai in (dir_db / 'fasta_by_segment').glob('*.fai'):
            data_length = pd.read_table(path_fai, usecols=[0, 1], names=['seq_id', 'length'])
            for seq_id, length in data_length.itertuples(index=False, name=None):
                len_by_segment[path_fai.name.split('.')[0]][seq_id] = length
        return len_by_segment

    def _merge_mash_output(self, mash_out_by_segment: dict[str, pd.DataFrame]) -> Path:
        """
        Merges the mash output for all segments.
        :param mash_out_by_segment: mash output by segment
        :return: None
        """
        for segment, data_mash in mash_out_by_segment.items():
            data_mash['segment'] = segment
        data_mash_combined = pd.concat([d for _, d in mash_out_by_segment.items()])
        path_out = self._folder / 'mash_combined.tsv'
        data_mash_combined.to_csv(path_out, sep='\t', index=False)
        logging.info(f'Combined mash output file created: {path_out.name} ({len(data_mash_combined):,} rows)')
        return path_out

    def _filter_mutually_exclusive_segments(
            self, ref_by_segment: dict[str, Union[pd.Series, None]], mut_exclusive_segments: list[str]) -> None:
        """
        Filters the mutually exclusive segments.
        :param ref_by_segment: Reference by segment
        :param mut_exclusive_segments: Mutually exclusive segments
        :return: Ref by segment
        """
        for keys_exclusive in [mut_exclusive_segments]:
            median_mult_by_key = [{
                'key': k,
                'dist': ref_by_segment[k].median_mult if ref_by_segment[k] is not None else 0
            } for k in keys_exclusive]
            median_mult_by_key.sort(key=lambda x: x['dist'], reverse=True)
            key_kept = median_mult_by_key[0]['key']
            keys_discarded = [x['key'] for x in median_mult_by_key[1:]]
            logging.info(f"Keeping '{key_kept}', discarding '{', '.join(keys_discarded)}'")
            for key in keys_discarded:
                ref_by_segment[key] = None

    def __create_fasta(self, ref_by_segment: dict[str, dict[str, Any]], dir_db: Path) -> Path:
        """
        Creates a FASTA file by combining the sequences of the selected references for each segment.
        :param ref_by_segment: Information on the selected reference genome by segment
        :param dir_db: Database directory
        :return: Path to merged FASTA
        """
        if not (dir_db / 'fasta_by_segment').exists():
            raise FileNotFoundError("Segment folder does not exist ('fasta_by_segment')")

        # Collect corresponding FASTA records
        records_out = []
        for segment, ref_info in ref_by_segment.items():
            if ref_info is None:
                continue
            path_fasta = dir_db / 'fasta_by_segment' / f'{segment}.fasta'
            command = Command(f"module load samtools; samtools faidx {path_fasta} {ref_info['ref_id']}")
            command.run(self.folder, disable_logging=True)
            if not command.returncode == 0:
                raise RuntimeError(f'Error extracting sequence: {command.stderr}')
            records_out.append(next(SeqIO.parse(StringIO(command.stdout), 'fasta')))

        # Save to merged FASTA file
        path_out = self.folder / 'merged_ref.fasta'
        with path_out.open('w') as handle:
            SeqIO.write(records_out, handle, 'fasta')
        logging.info(f'Merged FASTA file created ({len(records_out)} sequences)')
        return path_out
