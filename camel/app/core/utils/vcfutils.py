from pathlib import Path

import vcf

from camel.app.core.errors import InvalidParameterError
from camel.app.core.utils import fileutils


def is_multi_sample(vcf_file: Path) -> bool:
    """
    Function to check whether a VCF file contains multiple sample
    :param vcf_file: the vcf file to be checked (with complete path)
    :return: True if is a multiple sample VCF file
    """
    return len(vcf.Reader(filename=str(vcf_file)).samples) > 1

def retrieve_variants(vcf_file: Path, types: list[str] = None, excluded_types: list[str] = None):
    """
    Function to retrieve certain types of variants. Parameter
    'types' and 'excluded_types' is mutual exclusive. Either specify
    'types' to retrieved only certain variants or specify
    'excluded_types' to exclude variants.

    TYPES: 'unknown', ''indel', 'snp'
    SUBTYPES: 'tv', 'ts', 'ins', 'del'

    :param vcf_file: [optional] the vcf file to retrieve data
    :param types: [optional] types of variants to be retrieved
    :param excluded_types: types of variants to be excluded
    :return: list of records of types
    """
    types = [] if types is None else types
    excluded_types = [] if excluded_types is None else excluded_types

    if len(types) > 0 and len(excluded_types) > 0:
        raise InvalidParameterError(
            "Mutually exclusive parameters 'included types' and 'excluded types' are specified. Only one is "
            "allowed.")

    vcf_reader = vcf.Reader(filename=str(vcf_file))
    records = []
    if len(types) > 0:
        # only set types
        for rcd in vcf_reader:
            if any(x in types for x in [rcd.var_type, rcd.var_subtype]):
                records.append(rcd)
    elif len(excluded_types) > 0:
        # only set excluded_types
        for rcd in vcf_reader:
            if all(x not in excluded_types for x in [rcd.var_type, rcd.var_subtype]):
                records.append(rcd)
    else:
        # no conditions
        records = list(vcf_reader)
    return records

def count_variants(vcf_file: Path) -> int:
    """
    Counts the number of variants in a VCF file.
    :param vcf_file: VCF file.
    :return: Number of variants
    """
    if vcf_file is None:
        raise ValueError("VCF file should not be None")
    gzipped = fileutils.is_gzipped(vcf_file)
    with open(vcf_file, 'rb' if gzipped else 'r') as handle:
        vcf_reader = vcf.Reader(handle)
        variants = list(vcf_reader)
    return sum(variant.ALT != [None] for variant in variants if variant.FILTER is None or len(variant.FILTER) == 0)
