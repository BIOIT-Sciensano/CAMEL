import dataclasses
from pathlib import Path
from typing import Any

from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import vcfutils
from cyvcf2 import Variant

from camel.app.loggers import logger
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter


@dataclasses.dataclass
class FilterVariantsOutput:
    """
    Holder for the variant filtering output.
    """
    path_vcf: Path
    stats: dict[str, Any]
    informs: list[dict]


class FilterVariants:
    """
    Wrapper around the variant filters for the viral consensus pipeline.
    """

    PARAMS_BY_CALLER = {
        'bcftools': {
            'min_dp': {'expression': 'DP < {}', 'default': 10},
            'min_af': {'expression': '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])) < {}', 'default': 0.5},
            'min_qual': {'expression': 'QUAL < {}', 'default': 25},
            'min_mq': {'expression': 'MQ < {}', 'default': 30}
        },
        'clair3': {
            'min_dp': {'expression': 'DP < {}', 'default': 10},
            'min_af': {'expression': 'AF < {}', 'default': 0.5},
            'min_qual': {'expression': 'QUAL < {}', 'default': 25}
        }
    }

    def __init__(self, dir_: Path) -> None:
        """
        Initializes this workflow.
        :param dir_: Working directory
        :return: None
        """
        self._dir = dir_
        if not self._dir.exists():
            logger.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)

    def run(self, vcf_in: Path, calling_method: str, filters: dict[str, Any]) -> FilterVariantsOutput:
        """
        Runs the variant filtering workflow.
        :param vcf_in: Input VCF file
        :param calling_method: Variant calling method
        :param filters: Filters to apply, keyed by filter name
        :return: Variant filtering output
        """
        logger.info(f"Applying filters: {', '.join(filters.keys())}")
        path_vcf = vcf_in
        informs = []
        for filter_key, filter_value in filters.items():
            if filter_key not in FilterVariants.PARAMS_BY_CALLER[calling_method]:
                raise ValueError(f"Filter '{filter_key}' not supported for {calling_method}")
            logger.info(f'Applying filter: {filter_key} (value={filter_value})')
            path_vcf, inform = self._apply_filter(path_vcf, calling_method, filter_key, filter_value)
            informs.append(inform)
        return FilterVariantsOutput(path_vcf, FilterVariants._extract_stats(path_vcf), informs)

    def _apply_filter(self, path_vcf_in: Path, caller: str, filter_key: str, filter_value: Any) -> tuple[Path, dict]:
        """
        Applies the given variant filter to the input VCF file.
        :param path_vcf_in: Input VCF file
        :param caller: Variant caller
        :param filter_key: Filter key
        :param filter_value: Filter value
        :return: Path to filtered VCF file and tool informs
        """
        path_out = self._dir / f'variants-filt_{filter_key}.vcf'
        expression = str(FilterVariants.PARAMS_BY_CALLER[caller][filter_key]['expression']).format(filter_value)
        bcftools_filter = BcftoolsFilter()
        bcftools_filter.add_input_files({
            'VCF': [ToolIOFile(path_vcf_in)]
        })
        bcftools_filter.update_parameters(
            exclude=f'"{expression}"',
            output_filename=str(path_out),
            output_type='v',
            soft_filter=filter_key,
        )
        bcftools_filter.run()
        return path_out, bcftools_filter.informs

    @staticmethod
    def _extract_stats(path_vcf: Path) -> dict[str, Any]:
        """
        Extracts variant filtering stats by parsing the output VCF file.
        :param path_vcf: Input VCF file
        :return: Variant filtering statistics
        """
        variants: list[Variant] = vcfutils.parse_all_variants(path_vcf)
        return {
            'nb_variants': len(variants),
            'nb_snps': sum(v.is_snp for v in variants),
            'nb_snps_pass': sum(v.is_snp for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0)),
            'nb_indels': sum(v.is_indel for v in variants),
            'nb_indels_pass': sum(v.is_indel for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0)),
        }
