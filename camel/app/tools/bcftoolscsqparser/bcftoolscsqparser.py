import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as VcfRecord

from camel.app.toolkits.csq.csqutils import BCSQInfo
from camel.app.toolkits.csq.mutations.aminoacidmutation import AminoAcidMutation
from camel.app.toolkits.csq.mutations.basemutation import BaseMutation
from camel.app.toolkits.csq.mutations.frameshiftmutation import FrameshiftMutation
from camel.app.toolkits.csq.mutations.nucelotidemutation import NucleotideMutation
from camel.app.toolkits.csq.mutations.stopmutation import StopMutation
from camel.app.toolkits.csq.mutations.unknownmutation import UnknownMutation
from camel.app.toolkits.tabix import tabixparser
from camel.app.core.errors import InvalidToolInputError
from camel.app.core.io.tooliovalue import ToolIOValue
from camel.app.loggers import logger
from camel.app.core.tool import Tool


class CsqParsingError(ValueError):
    """
    Error that is raised when a mutation cannot be parsed.
    """
    pass


class CsqParser(Tool):
    """
    This tool is used to parse the annotated VCF files generated with bcftools csq.
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        """
        super().__init__('bcftools csq parser', '0.1')
        self._annot_tabix = {}

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF' not in self._tool_inputs:
            raise InvalidToolInputError("VCF input is required")
        if 'TSV' not in self._tool_inputs:
            logger.warning("TABIX annotation ('TSV') is missing, nucleotide mutations will be skipped.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        if 'TSV' in self._tool_inputs:
            self._annot_tabix = tabixparser.parse_tabix_annotation(self._tool_inputs['TSV'][0].path)

        with open(self._tool_inputs['VCF'][0].path) as handle:
            variants = list(vcf.Reader(handle))

        muts = []
        self._informs['unparsed'] = 0
        for v in variants:
            try:
                m = self.__parse_vcf_record(v)
                muts.append(m)
            except CsqParsingError as err:
                logger.warning(f"Unparsed mutation '{v}': {err}")
                self._informs['unparsed'] += 1

        self._informs['counts'] = {}
        for m in muts:
            key = m.__class__.__name__
            try:
                self._informs['counts'][key] += 1
            except KeyError:
                self._informs['counts'][key] = 1

        logger.info(f"Parsed CSQ mutations: {self._informs['counts']}")
        self._tool_outputs['VAL_mut'] = [ToolIOValue(m) for m in muts]

    def __parse_vcf_record(self, record: VcfRecord) -> BaseMutation:
        """
        Parses a mutation from a VCF record.
        :param record: Record
        :return: Parsed mutation
        """
        position = (record.CHROM, record.POS)
        if 'BCSQ' not in record.INFO:
            # When TABIX annotations are provided, the nucleotide mutations can be parsed
            if 'TSV' in self._tool_inputs:
                return NucleotideMutation.parse(record, self._annot_tabix)
            else:
                raise CsqParsingError(f'Mutation at position {position} does not contain BCSQ field: {record.INFO}')
        info = BCSQInfo.parse(record.INFO['BCSQ'][0])
        if info.type_ in ('missense', 'synonymous'):
            return AminoAcidMutation.parse(record, info)
        elif info.type_ == 'frameshift':
            return FrameshiftMutation.parse(record, info)
        elif info.type_ == 'non_coding':
            return NucleotideMutation.parse(record, self._annot_tabix)
        elif info.type_ == 'stop_gained':
            return StopMutation.parse(record, info)
        else:
            return UnknownMutation.parse(record, info)
