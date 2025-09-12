import argparse
from collections.abc import Sequence
from typing import Optional

import pandas as pd
from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.loggers import logger


class MainSnpDistances:
    """
    Script to determine SNP distances from a SNP matrix.
    """

    def __init__(self, args: Optional[Sequence[str]] = None):
        """
        Initializes the script.
        :param args: (Optional) arguments
        """
        self._args = MainSnpDistances._parse_args(args)

    @staticmethod
    def _parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments to parse
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', required=True, help='Input FASTA file')
        parser.add_argument('--output', required=True, help='Output file with SNP distances')
        return parser.parse_args(args)

    @staticmethod
    def calc_snp_distance(seq_a: str, seq_b: str) -> int:
        """
        Calculates the SNP distances between two sequences.
        N's are ignored in the calculation.
        :param seq_a: Sequence A
        :param seq_b: Sequence B
        :return: None
        """
        distance = 0
        for nuc_a, nuc_b in zip(seq_a, seq_b):
            if nuc_a == 'N' or nuc_b == 'N':
                continue
            distance += 1 if nuc_a != nuc_b else 0
        return distance

    def run(self) -> None:
        """
        Runs the script.
        :return: None
        """
        # Parse input
        with open(self._args.fasta) as handle:
            seqs = SeqIO.parse(handle, 'fasta')
            seq_by_id = {s.id: str(s.seq) for s in seqs}

        # Calculate distances
        data_distance = pd.DataFrame([{'sample': id_} for id_ in sorted(seq_by_id.keys())])
        for sample in sorted(seq_by_id.keys()):
            data_distance[sample] = data_distance.apply(
                lambda x: MainSnpDistances.calc_snp_distance(seq_by_id[x['sample']], seq_by_id[sample]), axis=1)

        # Save output
        data_distance.to_csv(self._args.output, sep='\t', index=False)
        logger.info(f"Output file created: {self._args.output}")


if __name__ == '__main__':
    Camel.get_instance()
    main = MainSnpDistances()
    main.run()
