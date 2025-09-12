import dataclasses
import logging
from pathlib import Path
from typing import Any

import vcf

from camel.app.command.command import Command


@dataclasses.dataclass
class FilterVariantsOutput:
    """
    Holder for the output of the variant filtering output.
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
            logging.info(f'Creating working directory: {self._dir}')
            self._dir.mkdir(parents=True)
        self._informs = []

    def run(self, vcf_in: Path, calling_method: str, filters: dict[str, Any]) -> FilterVariantsOutput:
        """
        Runs the variant filtering workflow.
        """
        logging.info(f"Applying filters: {', '.join(filters.keys())}")
        path_vcf = vcf_in
        for filter_key, filter_value in filters.items():
            if filter_key not in FilterVariants.PARAMS_BY_CALLER[calling_method]:
                raise ValueError(f"Filter '{filter_key}' not supported for {calling_method}")
            logging.info(f'Applying filter: {filter_key} (value={filter_value})')
            path_vcf = self.__apply_filter(path_vcf, calling_method, filter_key, filter_value)
        return FilterVariantsOutput(path_vcf, FilterVariants.__extract_stats(path_vcf), self._informs)

    def __apply_filter(self, path_vcf_in: Path, caller: str, filter_key: str, filter_value: str) -> Path:
        """
        Applies the given variant filter to the input VCF file.
        :param path_vcf_in: Input VCF file
        :param caller: Variant caller
        :param filter_key: Filter key
        :param filter_value: Filter value
        :return: Path to filtered VCF file
        """
        path_out = self._dir / f'variants-filt_{filter_key}.vcf'
        expression = str(FilterVariants.PARAMS_BY_CALLER[caller][filter_key]['expression']).format(filter_value)
        command = Command(' '.join([
            "module load bcftools/1.17;",
            f'bcftools filter --output-type v --soft-filter {filter_key} --exclude "{expression}" {path_vcf_in}',
            f'--output {path_out};'
        ]))
        command.run(self._dir)
        if not command.returncode == 0:
            raise RuntimeError(f'Error applying filter ({filter_key}): {command.stderr}')
        self._informs.append({'_name': 'bcftools filter 1.17', '_version': '1.17', '_command': command.command})
        return path_out

    @staticmethod
    def __extract_stats(path_vcf: Path) -> dict[str, Any]:
        """
        Extracts variant filtering stats by parsing the output VCF file.
        :param path_vcf: Input VCF file
        :return: Variant filtering statistics
        """
        variants = list(vcf.Reader(filename=str(path_vcf)))
        return {
            'nb_variants': len(variants),
            'nb_snps': sum(v.is_snp for v in variants),
            'nb_snps_filt': sum(v.is_snp for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0)),
            'nb_indels': sum(v.is_indel for v in variants),
            'nb_indels_filt': sum(v.is_indel for v in variants if (v.FILTER is None) or (len(v.FILTER) == 0)),
        }
