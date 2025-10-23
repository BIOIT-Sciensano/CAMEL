from dataclasses import dataclass

import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# noinspection PyProtectedMember
from vcf.model import _Record as VCFRecord

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


@dataclass(unsafe_hash=True, frozen=True)
class SNPPosition:
    """
    Represents a position in the SNP matrix.
    """
    contig: str
    position: int
    reference_base: str


class SnpMatrixConstructor(Tool):
    """
    Constructs SNP matrices in FASTA format based on VCF input files.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('SNP Matrix Constructor', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the tool input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError("No VCF input found.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        sample_names = [io.value for io in self._tool_inputs['SAMPLE_NAME']]
        logger.info(f"{len(sample_names)} samples provided")
        nucleotide_by_sample_by_position = self.__get_nucleotides_per_position()
        include_ref = 'include_ref' in self._parameters
        self._tool_outputs['FASTA'] = [self.__generate_matrix(
            sample_names, nucleotide_by_sample_by_position, include_ref)]

    def __get_nucleotides_per_position(self) -> dict[SNPPosition, dict[str, str]]:
        """
        Returns a dictionary with the nucleotide for each sample at each variant position.
        :return: Sample_names, nucleotide per position per sample
        """
        nucl_by_position = {}
        for vcf_file, io_sample in zip(self._tool_inputs['VCF'], self._tool_inputs['SAMPLE_NAME']):
            for record in SnpMatrixConstructor.parse_vcf_file(str(vcf_file.path), 'include_filtered' in self._parameters):
                position = SNPPosition(record.CHROM, record.POS, record.REF)
                if position not in nucl_by_position:
                    nucl_by_position[position] = {}
                if len(record.FILTER) > 0:
                    nucl_by_position[position][io_sample.value] = 'N'
                else:
                    nucl_by_position[position][io_sample.value] = str(record.ALT[0])

        # Remove positions with only N
        nucl_by_position = {pos: nucl for pos, nucl in nucl_by_position.items() if not all(
            [x == 'N' for x in nucl.values()])}

        logger.info(f"{len(nucl_by_position):,} SNP positions found across all samples")
        return nucl_by_position

    @staticmethod
    def parse_vcf_file(vcf_path: str, include_filtered: bool) -> list[VCFRecord]:
        """
        Parses a single VCF file.
        :param vcf_path: VCF path
        :param include_filtered: If True, filtered variants are included
        :return: List of SNP positions
        """
        vcf_records = []
        with open(vcf_path) as handle:
            for record in list(vcf.Reader(handle)):
                # Remove non SNP positions
                if not record.is_snp:
                    continue
                # Remove filtered records (unless specified otherwise)
                if (not include_filtered) and len(record.FILTER) > 0:
                    continue
                vcf_records.append(record)
        logger.info(f"{vcf_path} parsed: {len(vcf_records)} variant positions")
        return vcf_records

    def __generate_matrix(self, sample_names: list[str], nucl_by_pos: dict[SNPPosition, dict[str, str]],
                          include_ref: bool = True) -> ToolIOFile:
        """
        Generates a SNP matrix.
        :param sample_names: List of sample names
        :param nucl_by_pos: Nucleotide by sample by SNP position
        :param include_ref: If True, the reference is included in the SNP matrix
        :return: None
        """
        seq_by_sample_name = {name: [] for name in sample_names}
        if include_ref is True:
            seq_by_sample_name['reference'] = []
        for snp_pos, nucl_by_sample in sorted(nucl_by_pos.items(), key=lambda x: x[0].position):
            for name in seq_by_sample_name.keys():
                seq_by_sample_name[name].append(nucl_by_sample.get(name, snp_pos.reference_base))

        output_path = self.folder / 'snp_matrix.fasta'
        seqs = [SeqRecord(Seq(''.join(seq)), name, description='') for name, seq in seq_by_sample_name.items()]
        with open(output_path, 'w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        return ToolIOFile(output_path)
