from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from camelcore.app.io.tooliofile import ToolIOFile
from camelcore.app.utils import fastautils

from camel.app.core import toolutils
from camel.app.loggers import logger
from camel.app.tools.gatk import MASK_NT
from camel.app.tools.gatk.gatk import GATK


class GATKFastaAlternateReferenceMaker(GATK):
    """
    Generates an alternative reference sequence by applying the variants of a VCF file to the reference genome.

    Required inputs:
    - VCF: Variants to apply
    - FASTA_REF: Reference genome

    Optional inputs:
    - TXT_intervals: Intervals to extract (formatted as seq_id:start-end)
    - VCF_SNPmask: SNPs to mask

    Outputs:
    - FASTA: Alternative reference sequence(s)
    - FASTA_concatenated: Extracted intervals concatenated per reference sequence, with the gaps between the intervals
      masked (only if 'concatenate_sequence_segments' is set)
    """

    def __init__(self) -> None:
        """
        Initializes this tool.
        :return: None
        """
        super().__init__('gatk FastaAlternateReferenceMaker')

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        toolutils.check_input(self, keys_required=['VCF', 'FASTA_REF'])
        if self.get_param_value('concatenate_sequence_segments') and 'TXT_intervals' not in self._tool_inputs:
            logger.warning("Concatenating sequence segments requires the 'TXT_intervals' input, option disabled")
            self.update_parameters(concatenate_sequence_segments=False)
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        super()._execute_tool()
        if self.get_param_value('concatenate_sequence_segments'):
            self.__concatenate_sequence_segments()

    def _get_input_options(self) -> list[str]:
        """
        Returns the command line options for the input files.
        :return: Input options
        """
        options = [*super()._get_input_options(), f"-V {self._tool_inputs['VCF'][0].path}"]
        if 'VCF_SNPmask' in self._tool_inputs:
            options.append(f"--snpmask {self._tool_inputs['VCF_SNPmask'][0].path}")
        return options

    def _get_parameter_options(self) -> list[str]:
        """
        Returns the command line options for the parameters.
        :return: Parameter options
        """
        return self._build_options(excluded_parameters=['concatenate_sequence_segments'])

    def _set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self.folder / self.get_param_value('output'))]
        if self.get_param_value('concatenate_sequence_segments'):
            self._tool_outputs['FASTA_concatenated'] = [
                ToolIOFile(self.folder / 'consensus_sequences_concatenated.fa')]
        else:
            self._tool_outputs['FASTA_concatenated'] = []

    def __read_intervals(self) -> list[tuple[str, int, int]]:
        """
        Reads the intervals and sorts them the same way as GATK does (by reference sequence and start position).
        :return: List of intervals (sequence id, start, end)
        """
        with self._tool_inputs['FASTA_REF'][0].path.open() as handle:
            seq_order = {record.id: idx for idx, record in enumerate(SeqIO.parse(handle, 'fasta'))}
        logger.debug(f'Reference sequence ids: {list(seq_order.keys())}')

        intervals = []
        with self._tool_inputs['TXT_intervals'][0].path.open() as handle:
            for line in handle:
                if line.strip() == '':
                    continue
                seq_id, positions = line.strip().rsplit(':', 1)
                start, end = positions.split('-')
                intervals.append((seq_id, int(start), int(end)))
        return sorted(intervals, key=lambda interval: (seq_order[interval[0]], interval[1]))

    def __concatenate_sequence_segments(self) -> None:
        """
        Concatenates the extracted sequence segments per reference sequence, masking the gaps between the segments.
        The extracted sequences are relabeled with their interval.
        :return: None
        """
        path_extracted = self._tool_outputs['FASTA'][0].path
        records_extracted = fastautils.read_as_dict(path_extracted)

        # GATK labels the extracted sequences with their (1-based) index in the sorted intervals
        segments = {}
        last_end = {}
        for idx, (seq_id, start, end) in enumerate(self.__read_intervals(), start=1):
            record = records_extracted[str(idx)]
            record.id = f'{seq_id}:{start}-{end}'
            if seq_id in segments:
                segments[seq_id].append(MASK_NT * (start - last_end[seq_id] - 1))
            else:
                segments[seq_id] = []
            segments[seq_id].append(str(record.seq))
            last_end[seq_id] = end

        records_concatenated = [
            SeqRecord(Seq(''.join(parts)), id=seq_id, description=seq_id) for seq_id, parts in segments.items()]
        logger.debug(f'Concatenated sequences: {[f"{x.id}/{len(x.seq)}" for x in records_concatenated]}')
        fastautils.write(list(records_extracted.values()), path_extracted)
        fastautils.write(records_concatenated, self._tool_outputs['FASTA_concatenated'][0].path)
