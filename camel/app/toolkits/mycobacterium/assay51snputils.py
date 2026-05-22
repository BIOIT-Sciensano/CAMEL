from dataclasses import dataclass
from pathlib import Path

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


@dataclass(frozen=True, unsafe_hash=True)
class SNPPosition:
    """
    Holder class for SNP positions.
    """
    pos: int
    name: str
    ref: str

    alt_unfilt: str | None = None
    alt_filt: str | None = None
    is_unfilt_snp: bool = False
    is_filt_snp: bool = False

    @property
    def nucl(self) -> str:
        """
        Returns the nucleotide at the given position.
        :return: Nucleotide at the given position
        """
        if self.is_filt_snp and self.alt_filt:
            return self.alt_filt
        elif self.is_unfilt_snp and self.alt_unfilt:
            return self.alt_unfilt
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
        if not self.is_filt_snp and not self.is_unfilt_snp:
            return 'green'
        elif self.is_filt_snp and self.is_unfilt_snp:
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
            positions.append(SNPPosition(pos=int(parts[2]), name=parts[3], ref=parts[4]))
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
