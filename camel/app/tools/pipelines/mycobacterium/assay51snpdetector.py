import dataclasses

from camelcore.app.utils import vcfutils

from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.app.toolkits.mycobacterium import assay51snputils
from camel.app.toolkits.mycobacterium.assay51snputils import SCGProfile, SNPPosition


class Assay51SnpDetector(Tool):
    """
    This class reports the detected species based on a set of 51 SNPs.
    The SNPs and the methodology are described in:
    https://jcm.asm.org/content/52/6/1962.full
    """

    SCG_SNPS_RANGE = range(7, 52)

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Mycobacterium: 51SNP detector', '0.1')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Parse input files
        snp_positions = assay51snputils.parse_snp_positions(
            self._tool_inputs['BED'][0].path
        )
        snp_positions = self.__add_vcf_info(snp_positions)
        snp_position_by_name = {p.name: p for p in snp_positions}
        scg_profiles = assay51snputils.parse_scg_profiles(
            self._tool_inputs['TSV'][0].path
        )

        # Extract required info
        self.__extract_positive_control_snp(snp_position_by_name)
        self.__extract_gyrb_group(snp_position_by_name)
        self.__extract_genetic_group(snp_position_by_name)
        self.__extract_best_scg_group(snp_position_by_name, scg_profiles)
        self.__extract_snp_informs(snp_position_by_name)

        # Update the informs
        self._informs['snp_positions_by_name'] = {
            name: dataclasses.asdict(snp) for name, snp in snp_position_by_name.items()
        }

    def __parse_vcf_variants(self, vcf_key: str) -> dict:
        """
        Parses VCF variants and returns a dict indexed by position.
        All variants (including indels/MNPs) are retained; only the first nucleotide of the ALT allele is
        used downstream (see __get_alt_allele), which keeps the single-base, position-wise comparison the
        51 SNP assay relies on while staying backwards-equivalent with the original behaviour.
        :param vcf_key: Key for tool inputs ('VCF' or 'VCF_filt')
        :return: Dictionary mapping positions to variant records
        """
        variants = vcfutils.parse_all_variants(self._tool_inputs[vcf_key][0].path)
        records_by_pos = {v.POS: v for v in variants if len(v.ALT) > 0}
        logger.debug(f'{len(records_by_pos)} variants parsed from {vcf_key}')
        return records_by_pos

    @staticmethod
    def __get_alt_allele(records_by_pos: dict, pos: int) -> str | None:
        """
        Safely extracts the ALT allele from variant records.
        Only the first nucleotide is returned, so indels/MNPs collapse to a single base. This matches the
        original behaviour of the 51 SNP assay (backwards equivalence).
        :param records_by_pos: Dictionary of variant records by position
        :param pos: Position to look up
        :return: First nucleotide of the ALT allele as string or None
        """
        return str(records_by_pos[pos].ALT[0])[0] if pos in records_by_pos else None

    def __update_snp_position(self, snp_position: SNPPosition, records_by_pos: dict,
                              records_filt_by_pos: dict) -> SNPPosition:
        """
        Updates a SNP position with VCF information.
        :param snp_position: Original SNP position
        :param records_by_pos: All variants indexed by position
        :param records_filt_by_pos: Filtered variants indexed by position
        :return: Updated SNP position
        """
        return dataclasses.replace(
            snp_position,
            alt_unfilt=self.__get_alt_allele(records_by_pos, snp_position.pos),
            is_unfilt_snp=snp_position.pos in records_by_pos,
            alt_filt=self.__get_alt_allele(records_filt_by_pos, snp_position.pos),
            is_filt_snp=snp_position.pos in records_filt_by_pos,
        )

    def __add_vcf_info(self, snp_positions: list[SNPPosition]) -> list[SNPPosition]:
        """
        Adds VCF records to the SNP positions, if no SNP is found the value is kept at None.
        :param snp_positions: SNP positions
        :return: List of updated SNP positions
        """
        records_by_pos = self.__parse_vcf_variants('VCF')
        records_filt_by_pos = self.__parse_vcf_variants('VCF_filt')
        return [
            self.__update_snp_position(snp_position, records_by_pos, records_filt_by_pos)
            for snp_position in snp_positions
        ]

    def __extract_positive_control_snp(
        self, snp_positions_by_name: dict[str, SNPPosition]
    ) -> None:
        """
        Checks the positive control SNP position.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['mtbc_pos_control'] = (
            'OK'
            if (
                not snp_positions_by_name['SNP01'].is_filt_snp
                and not snp_positions_by_name['SNP01'].is_unfilt_snp
            )
            else 'NOT OK'
        )

    def __extract_gyrb_group(
        self, snp_positions_by_name: dict[str, SNPPosition]
    ) -> None:
        """
        Extracts the gyrB group.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['gyrB_group'] = 'NA'
        self._informs['gyrB_species'] = 'NA'
        for profile in assay51snputils.GYRB_PROFILES:
            if all(
                snp_positions_by_name[key].nucl == profile[key]
                for key in ['SNP02', 'SNP03', 'SNP04']
            ):
                self._informs['gyrB_group'] = profile['group']
                self._informs['gyrB_species'] = profile['species']

    def __extract_genetic_group(
        self, snp_positions_by_name: dict[str, SNPPosition]
    ) -> None:
        """
        Extracts the genetic group.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['genetic_group'] = 'NA'
        for genetic_group in assay51snputils.GENETIC_GROUPS:
            if all(
                snp_positions_by_name[key].nucl == genetic_group[key]
                for key in ['SNP05', 'SNP06']
            ):
                self._informs['genetic_group'] = genetic_group['name']

    def __extract_best_scg_group(
        self, snp_positions_by_name: dict[str, SNPPosition], profiles: list[SCGProfile]
    ) -> None:
        """
        Returns the best matching SCG profile.
        :param snp_positions_by_name: SNP positions by name
        :param profiles: List of SCGprofiles
        :return: Best matching profile
        """
        nucl_str = ''.join(
            snp_positions_by_name[f'SNP{i:02d}'].nucl
            for i in self.SCG_SNPS_RANGE
        )
        best_profile = max(
            profiles,
            key=lambda p: sum(p.snps[i] == nucl_str[i] for i in range(len(nucl_str)))
        )
        nb_matches = sum(best_profile.snps[i] == nucl_str[i] for i in range(len(nucl_str)))
        self._informs['scg_profile'] = dataclasses.asdict(best_profile)
        self._informs['scg_nb_snps_matched'] = nb_matches

    def __extract_snp_informs(
        self, snp_positions_by_name: dict[str, SNPPosition]
    ) -> None:
        """
        Extracts the informs for the individual SNPs.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        for name, pos in snp_positions_by_name.items():
            if pos.is_filt_snp:
                self._informs[name] = pos.alt_filt
            elif pos.is_unfilt_snp:
                self._informs[name] = f'{pos.alt_unfilt}*'
            else:
                self._informs[name] = '-'
