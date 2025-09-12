from dataclasses import dataclass
from pathlib import Path
# noinspection PyProtectedMember
from vcf.model import _Record as Record

GYRB_PROFILES = [
    {'group': 'TC1', 'SNP02': 'G', 'SNP03': 'G', 'SNP04': 'C', 'species':
        '<i>M. tuberculosis, M. africanum II, M. canettii</i>'},
    {'group': 'TC2', 'SNP02': 'T', 'SNP03': 'A', 'SNP04': 'C', 'species':
        '<i>M. bovis, M.bovis BCG, M. caprae</i>'},
    {'group': 'TC3', 'SNP02': 'T', 'SNP03': 'G', 'SNP04': 'C', 'species': '<i>M. piniipedii, M. africanum I</i>'},
    {'group': 'TC4', 'SNP02': 'T', 'SNP03': 'G', 'SNP04': 'T', 'species': '<i>M. microti</i>'}
]

GENETIC_GROUPS = [
    {'name': 'GG1', 'SNP05': 'T', 'codon05': 'CTG (Leu)', 'SNP06': 'C', 'codon06': 'ACC (Thr)'},
    {'name': 'GG2', 'SNP05': 'G', 'codon05': 'CGG (Arg)', 'SNP06': 'C', 'codon06': 'ACC (Thr)'},
    {'name': 'GG3', 'SNP05': 'G', 'codon05': 'CGG (Arg)', 'SNP06': 'G', 'codon06': 'AGC (Ser)'}
]


@dataclass
class SNPPosition:
    """
    Holder class for SNP positions.
    """
    pos: int
    name: str
    ref: str
    vcf_record: Record = None
    vcf_filt_record: Record = None

    @property
    def nucl(self) -> str:
        """
        Returns the nucleotide at the given position.
        :return: Nucleotide at the given position
        """
        if (self.vcf_filt_record is not None) and self.vcf_filt_record.is_snp:
            return str(self.vcf_filt_record.ALT[0])
        elif (self.vcf_record is not None) and self.vcf_record.is_snp:
            return str(self.vcf_record.ALT[0])
        else:
            return self.ref

    def get_color(self, ref: str) -> str:
        """
        Returns the color for the SNP cell.
        If the position does not match the reference the returned color is 'red'. If the position matches based
        on the filtered VCF file the returned color is 'green'. If only the unfiltered VCF supports the nucleotide
        'lightgreen' is returned.
        :param ref: Reference base
        :return: Color
        """
        if self.nucl != ref:
            return 'red'
        elif (self.vcf_filt_record is None) and (self.vcf_record is None):
            return 'green'
        elif (self.vcf_filt_record is not None) and (self.vcf_record is not None):
            return 'green'
        else:
            return 'lightgreen'


@dataclass
class SCGProfile:
    """
    Holder class for SCG profiles.
    """
    st: str
    scg: str
    snps: str


def parse_snp_positions(positions_path: Path) -> list[SNPPosition]:
    """
    Parses the SNP positions from the tabular input file.
    :param positions_path: Path to the SNP positions BED file
    :return: List of parsed positions
    """
    positions = []
    with open(positions_path) as handle:
        for line in handle.readlines():
            parts = line.strip().split('\t')
            positions.append(SNPPosition(int(parts[2]), parts[3], parts[4]))
    return sorted(positions, key=lambda p: p.name)


def parse_scg_profiles(profiles_path: Path) -> list[SCGProfile]:
    """
    Parses the SNP cluster group (SCG) profiles.
    :param profiles_path: Profiles file path
    :return: List of profiles
    """
    profiles = []
    with open(profiles_path) as handle:
        for line in handle.readlines()[1:]:
            parts = line.strip().split('\t')
            profiles.append(SCGProfile(parts[0], parts[1], ''.join(parts[2:])))
    return profiles
