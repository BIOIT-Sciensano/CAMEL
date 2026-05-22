import hashlib
import json
import logging
from pathlib import Path

from Bio import SeqIO
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool


class CollectIterativeMappingStats(Tool):
    """
    Collects stats for the iterative mapping.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Collect iterative mapping stats', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidToolInputError("Consensus FASTA input is required ('FASTA')")
        if 'JSON_depth' not in self._tool_inputs:
            raise InvalidToolInputError("Depth information is required ('JSON_depth')")
        if 'VCF_p1' not in self._tool_inputs:
            raise InvalidToolInputError("Phase 1 VCF input is required ('VCF_p1')")
        if 'VCF_p2' not in self._tool_inputs:
            raise InvalidToolInputError("Phase 2 VCF input is required ('VCF_p2')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse depth information
        with self._tool_inputs['JSON_depth'][0].path.open() as handle:
            data_depth = json.load(handle)
        seq_ids = list(data_depth['by_chr'].keys())
        self._informs['seq_ids'] = seq_ids

        # Initialize informs
        self._informs['all_segments'] = {
            'iter': self._parameters['nb_iter'].value,
            'dirname': next(p.name for p in self._tool_inputs['FASTA'][0].path.parents if p.name.startswith('iter_'))
        }
        self._informs['by_segment'] = {seq_id: {} for seq_id in seq_ids}

        # Add consensus sequence hash
        self._informs['all_segments']['sequence_md5'] = CollectIterativeMappingStats._hash_consensus_sequence(
            self._tool_inputs['FASTA'][0].path)

        # Add depth information
        for k, v in data_depth['total'].items():
            self._informs['all_segments'][k] = v
        for seq_id, stat_by_key in data_depth['by_chr'].items():
            for k, v in stat_by_key.items():
                self._informs['by_segment'][seq_id][k] = v

        # Parse VCF files
        self._collect_variant_stats(self._tool_inputs['VCF_p1'][0].path, 1, seq_ids)
        self._collect_variant_stats(self._tool_inputs['VCF_p2'][0].path, 2, seq_ids)

        # Export to JSON
        self._export_informs_to_json()

    @staticmethod
    def _hash_consensus_sequence(path_fasta: Path) -> str:
        """
        Hashes the current consensus sequence to track if it changed.
        :param path_fasta: Consensus sequence FASTA file
        :return: Hash
        """
        # Calculate the hash of the sequence
        with path_fasta.open() as handle:
            full_seq = '__'.join(str(s.seq) for s in SeqIO.parse(handle, 'fasta'))
            return hashlib.md5(full_seq.encode('ascii')).hexdigest()

    @staticmethod
    def _count_variant_types(variants: list) -> dict[str, int]:
        """
        Counts the different variant types in the input list.
        :param variants: List of variants
        :return: Dictionary with nb. of variants by type
        """
        return {
            'nb_variants': len(variants),
            'nb_snps': sum(v.is_snp for v in variants),
            'nb_snps_filt': sum(v.is_snp for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0)),
            'nb_indels': sum(v.is_indel for v in variants),
            'nb_indels_filt': sum(v.is_indel for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0))
        }

    def _collect_variant_stats(self, path_vcf: Path, phase: int, seq_ids: list[str]) -> None:
        """
        Collects stats concerning the detected variants.
        :param path_vcf: Input VCF file
        :param phase: Current phase
        :param seq_ids: List of sequence ids
        :return: None
        """
        # Parse VCF input
        variants = vcfutils.parse_all_variants(path_vcf)

        # Group variants by seq id
        variants_by_seq_id = {seq_id: [] for seq_id in seq_ids}
        for v in variants:
            if v.CHROM not in variants_by_seq_id:
                variants_by_seq_id[v.CHROM] = []
            variants_by_seq_id[v.CHROM].append(v)

        # Collect informs
        for k, v in CollectIterativeMappingStats._count_variant_types(variants).items():
            self._informs['all_segments'][f'phase_{phase}-{k}'] = v

        # By segment
        for seq_id in seq_ids:
            for k, v in CollectIterativeMappingStats._count_variant_types(variants_by_seq_id[seq_id]).items():
                self.informs['by_segment'][seq_id][f'phase-{phase}-{k}'] = v

    def _export_informs_to_json(self) -> None:
        """
        Exports the informs to JSON format.
        """
        path_out = self.folder / self._parameters['output_filename'].value
        with path_out.open('w') as handle:
            json.dump(self.informs, handle, indent=2)
        logging.info(f'Statistics exported to: {path_out}')
        self._tool_outputs['JSON'] = [ToolIOFile(path_out)]
