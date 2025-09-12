from typing import Any

import pandas as pd
import vcf

from camel.app.error import InvalidToolInputError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class CallMultiAllelicSites(Tool):
    """
    Calls multi-allelic sites in the input pileup in VCF format.
    """

    AMBIGUITY_BASES = {'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K'}

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('Call multi-allelic sites', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError('Pileup input is required (VCF)')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        records_out = []
        logger.info(f"Calling multi-allelic sites from: {self._tool_inputs['VCF'][0].path}")
        with self._tool_inputs['VCF'][0].path.open() as handle:
            for variant in vcf.Reader(handle):
                # Skip invariant & low depth sites
                if len(variant.ALT) == 1 or variant.INFO['DP'] < int(self._parameters['min_dp'].value):
                    continue

                # Calculate allele frequencies
                allele_freq_nucleotide = [(str(allele), dp / sum(variant.samples[0].data.AD)) for allele, dp in zip(
                    variant.alleles, variant.samples[0].data.AD)]
                allele_freq_nucleotide.sort(key=lambda x: x[-1], reverse=True)

                # Get the majority allele
                allele_major, freq = allele_freq_nucleotide[0]
                if allele_major not in 'ACTG':
                    raise ValueError(f'Invalid majority allele: {allele_major}')

                # Check the frequency of the most common alternate allele
                minor_allele_base, minor_allele_freq = allele_freq_nucleotide[1]
                if minor_allele_freq < float(self._parameters['min_freq_minor_allele'].value):
                    continue

                key = ''.join(sorted([allele_major, allele_freq_nucleotide[1][0]]))
                records_out.append({
                    'chrom': variant.CHROM,
                    'pos': variant.POS,
                    'major_allele': allele_major,
                    'secondary_allele': allele_freq_nucleotide[1][0] if allele_freq_nucleotide[1][1] > 0.33 else 'N',
                    'major_freq': freq,
                    'dp': variant.INFO['DP'],
                    'alleles': ','.join(str(x) for x in variant.alleles),
                    'ad': ','.join(str(x) for x in variant.samples[0].data.AD),
                    'iupac': CallMultiAllelicSites.AMBIGUITY_BASES[key]
                })
        self._set_output(records_out)

    def _set_output(self, records_out: list[dict[str, Any]]) -> None:
        """
        Collects the tool output.
        :param records_out: Output records with multi-allelic sites.
        """
        path_out = self.folder / 'multi_allelic_sites.tsv'
        pd.DataFrame(records_out).to_csv(path_out, sep='\t', index=False)
        self._tool_outputs['TSV'] = [ToolIOFile(path_out)]
        self._informs['nb_sites'] = len(records_out)
        logger.info(f'{len(records_out):,} multi-allelic sites detected')
        self._informs['min_freq_minor_allele'] = self._parameters['min_freq_minor_allele'].value
        self._informs['min_dp'] = int(self._parameters['min_dp'].value)
