from pathlib import Path
from typing import Dict, Any

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

OUTPUT_STATS = 'variant_filtering/stats/json.io'
OUTPUT_INFORMS_ALL = 'variant_filtering/informs_all.io'
OUTPUT_SUMMARY = 'variant_filtering/summary.{ext}'
OUTPUT_VCF = 'variant_filtering/06-regions/vcf.io'


def get_filtering_param(config: dict[str, Any], filter_key: str, param_name: str) -> Any:
    """
    Returns the value of the given filtering parameter, returns 'None' when the parameter is not set.
    :param config: Snakemake configuration
    :param filter_key: Filter key
    :param param_name: Parameter name
    :return: Parameter value if set, None otherwise
    """
    if filter_key not in config['variant_filtering']:
        return None
    return config['variant_filtering'][filter_key].get(param_name)


def filter_is_disabled(filter_key: str, config: Dict[str, Any]) -> bool:
    """
    Returns True if the given filter is disabled, False otherwise.
    :param filter_key: Filter key
    :param config: Snakemake config
    :return: True if disabled, false otherwise
    """
    if filter_key not in config['variant_filtering']:
        return False
    return config['variant_filtering'][filter_key].get('disabled', False)
