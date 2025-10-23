
from vcf import VCFReader

from camel.app.toolkits.mycobacterium import assay51snputils
from camel.app.toolkits.mycobacterium.assay51snputils import SNPPosition, SCGProfile
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class Assay51SnpDetector(Tool):
    """
    This class reports the detected species based on a set of 51 SNPs.
    The SNPs and the methodology are described in:
    https://jcm.asm.org/content/52/6/1962.full
    """

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
        snp_positions = assay51snputils.parse_snp_positions(self._tool_inputs['BED'][0].path)
        self._informs['snp_positions_by_name'] = {p.name: p for p in snp_positions}
        scg_profiles = assay51snputils.parse_scg_profiles(self._tool_inputs['TSV'][0].path)
        self.__add_vcf_records(snp_positions)

        # Extract required info
        self.__extract_positive_control_snp(self._informs['snp_positions_by_name'])
        self.__extract_gyrb_group(self._informs['snp_positions_by_name'])
        self.__extract_genetic_group(self._informs['snp_positions_by_name'])
        self.__extract_best_scg_group(self._informs['snp_positions_by_name'], scg_profiles)
        self.__extract_snp_informs(self._informs['snp_positions_by_name'])

    def __add_vcf_records(self, snp_positions: list[SNPPosition]) -> None:
        """
        Adds VCF records to the SNP positions, if no SNP is found the value is kept at None.
        :param snp_positions: SNP positions
        :return: None
        """
        with open(self._tool_inputs['VCF'][0].path) as handle:
            records_by_pos = {v.POS: v for v in VCFReader(handle)}
            logger.debug(f'{len(records_by_pos)} SNPs parsed')
        with open(self._tool_inputs['VCF_filt'][0].path) as handle:
            records_filt_by_pos = {v.POS: v for v in VCFReader(handle)}
            logger.debug(f'{len(records_filt_by_pos)} filtered SNPs parsed')
        for snp_position in snp_positions:
            if snp_position.pos in records_by_pos:
                snp_position.vcf_record = records_by_pos[snp_position.pos]
            if snp_position.pos in records_filt_by_pos:
                snp_position.vcf_filt_record = records_filt_by_pos[snp_position.pos]

    def __extract_positive_control_snp(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Checks the positive control SNP position.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['mtbc_pos_control'] = 'OK' if (
                snp_positions_by_name['SNP01'].vcf_record is None and
                snp_positions_by_name['SNP01'].vcf_filt_record is None) else 'NOT OK'

    def __extract_gyrb_group(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Extracts the gyrB group.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['gyrB_group'] = 'NA'
        self._informs['gyrB_species'] = 'NA'
        for profile in assay51snputils.GYRB_PROFILES:
            if all(snp_positions_by_name[key].nucl == profile[key] for key in ['SNP02', 'SNP03', 'SNP04']):
                self._informs['gyrB_group'] = profile['group']
                self._informs['gyrB_species'] = profile['species']

    def __extract_genetic_group(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Extracts the genetic group.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        self._informs['genetic_group'] = 'NA'
        for genetic_group in assay51snputils.GENETIC_GROUPS:
            if all(snp_positions_by_name[key].nucl == genetic_group[key] for key in ['SNP05', 'SNP06']):
                self._informs['genetic_group'] = genetic_group['name']

    def __extract_best_scg_group(self, snp_positions_by_name: dict[str, SNPPosition], profiles: list[SCGProfile]) ->\
            None:
        """
        Returns the best matching SCG profile.
        :param snp_positions_by_name: SNP positions by name
        :param profiles: List of SCGprofiles
        :return: Best matching profile
        """
        keys = ['SNP{:02d}'.format(i) for i in range(7, 52)]
        nucl_str = ''.join([snp_positions_by_name[k].nucl for k in keys])
        counts = []
        for profile in profiles:
            nb_matching_snps = sum([profile.snps[i] == nucl_str[i] for i in range(0, len(nucl_str))])
            counts.append((profile, nb_matching_snps))
        self._informs['scg_profile'], self._informs['scg_nb_snps_matched'] = sorted(counts, key=lambda x: -x[-1])[0]

    def __extract_snp_informs(self, snp_positions_by_name: dict[str, SNPPosition]) -> None:
        """
        Extracts the informs for the individual SNPs.
        :param snp_positions_by_name: SNP positions by name
        :return: None
        """
        for name, pos in snp_positions_by_name.items():
            if (pos.vcf_record is None) and (pos.vcf_filt_record is None):
                self._informs[name] = '-'
            elif pos.vcf_filt_record is None:
                self._informs[name] = '{}*'.format(str(pos.vcf_record.ALT[0]))
            else:
                self._informs[name] = str(pos.vcf_filt_record.ALT[0])
