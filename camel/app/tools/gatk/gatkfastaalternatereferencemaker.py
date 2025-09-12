from Bio import SeqIO, Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from camel.app.components.files.fastautils import FastaUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.gatk import MASK_NT
from camel.app.tools.gatk.gatk import GATK


class GATKFastaAlternateReferenceMaker(GATK):
    """
    Class for GATK FastaAlternateReferenceMaker function
    """

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('gatk FastaAlternateReferenceMaker', '3.7')
        self._required_inputs = ['VCF']
        self._output_type = 'FASTA'
        self._specific_parameters = ['concatenate_sequence_segments']
        self._fasta_concatenated = 'consensus_sequences_concatenated.fa'
        self._fasta_extracted = None
        self._concatenate_sequence = False

    def _execute_tool(self):
        """
        Run GATK FastaAlternateReferenceMaker
        :return: None
        """
        super(GATKFastaAlternateReferenceMaker, self)._execute_tool()

        if self._concatenate_sequence:
            FastaUtils.write(self.__concatenate_sequence_segments(), self._fasta_concatenated)

    def _check_parameters(self):
        """
        Checks tool parameters
        :return: None
        """
        super(GATKFastaAlternateReferenceMaker, self)._check_parameters()

        if 'concatenate_sequence_segments' in self._parameters:
            self._concatenate_sequence = True

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        if self._concatenate_sequence:
            if 'TXT_intervals' not in self._tool_inputs:
                logger.warning(
                    "FastaAlternateReferenceMaker opt 'concatenate_sequence_segments' required 'TXT_intervals' input is missing, option disabled.")
                self._concatenate_sequence = False

        super(GATKFastaAlternateReferenceMaker, self)._check_input()

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKFastaAlternateReferenceMaker, self)._set_input()

        if 'VCF_SNPmask' in self._tool_inputs:
            self._input_string += "--snpmask {} ".format(self._tool_inputs['VCF_SNPmask'][0].path)

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        super(GATKFastaAlternateReferenceMaker, self)._set_output()
        # set default output type self._output_type: 'FASTA'
        self._fasta_extracted = self._tool_outputs['FASTA'][0].path

        if self._concatenate_sequence:
            # when 'concatenate_sequence_segments' set, final FASTA output will be
            # generated from FastaAlternateReferenceMaker output
            self._fasta_concatenated = self._folder / self._fasta_concatenated
            self._tool_outputs['FASTA_concatenated'] = [ToolIOFile(self._fasta_concatenated)]
        else:
            self._tool_outputs['FASTA_concatenated'] = []

    @staticmethod
    def __get_interval_inform(interval):
        """
        Extract interval information from a interval specificaiton seqid:pstart-pend
        :return: seqid, sequence id
        :return: pstart, start position
        :return: pend, end position
        """
        seq_id, interval_pos = interval.split(':')
        pstart, pend = interval_pos.split('-')

        return seq_id, int(pstart), int(pend)

    def __rearrange_seq_intervals(self, seq_intervals):
        """
        Rearrange seq intervals according to the FASTA_REF seqid order, as GATK keeps sequence in this ordering
        :return: seq intervals ordered (list)
        """
        seq_grouped_intervals = {}
        for idx, interval in enumerate(seq_intervals):
            seq_id, pstart, pend = self.__get_interval_inform(interval)
            if seq_id in seq_grouped_intervals:
                seq_grouped_intervals[seq_id].append(interval)
            else:
                seq_grouped_intervals[seq_id] = [interval]

        refseq_ids = [x.id for x in list(SeqIO.parse(self._tool_inputs['FASTA_REF'][0], "fasta"))]

        logger.debug("  refseq_ids: {}".format(refseq_ids))

        seq_intervals_ordered = []
        for seqid in refseq_ids:
            if seqid in seq_grouped_intervals:
                seq_intervals_ordered += seq_grouped_intervals[seqid]

        return seq_intervals_ordered

    def __concatenate_sequence_segments(self):
        """
        Concatenate sequence segments for each sequence (identified by seqid from interval file)
        :return: segments concatenated sequences
        """
        extracted_seq_dict = FastaUtils.read_as_dict(self._fasta_extracted)
        seq_intervals = [x.strip() for x in open(self._tool_inputs['TXT_intervals'][0].path, 'r')]
        seq_intervals_ordered = self.__rearrange_seq_intervals(seq_intervals)

        concatenated_seqs = []
        last_seq_id = None
        last_end_pos = 0
        concatenate_seq = ''
        # Note:
        # - the extracted sequences are based on seq_intervals, hence there is 1-to-1 map
        # - GATK FastaAlternateReferenceMaker generate sequence with number as id (1-based)
        for idx, interval in enumerate(seq_intervals_ordered):
            seq_record = extracted_seq_dict[str(idx + 1)]
            seq_record.id = interval
            seq_id, pos_start, pos_end = self.__get_interval_inform(interval)
            if seq_id == last_seq_id:
                # from same sequence segament
                # add (mask sequence (Ns) + new segment)
                numb_padding = pos_start - last_end_pos - 1
                logger.debug("number of padding: {}, seq_length {}, interval {}, seq_record {}".format(
                    numb_padding, len(seq_record.seq), interval, str(seq_record)))
                # noinspection PyTypeChecker
                concatenate_seq += MASK_NT * numb_padding + str(seq_record.seq)
            else:
                # a new sequence segament
                if last_seq_id is not None:
                    concatenated_seqs.append(
                        SeqRecord(Seq(concatenate_seq, Alphabet), id=last_seq_id,
                                  description=last_seq_id))
                concatenate_seq = str(seq_record.seq)
                last_seq_id = seq_id
            last_end_pos = pos_end

        # add the last sequence
        concatenated_seqs.append(SeqRecord(Seq(concatenate_seq, Alphabet), id=last_seq_id, description=last_seq_id))
        logger.debug(" concatenated seqs inform: {}".format(
            ["{}/{}".format(x.id, len(x.seq)) for x in concatenated_seqs]))

        # Note: seq_record in extracted_seq_dict has been updated with ids, now output into self._fasta_extracted with
        #       updated ids
        FastaUtils.write(list(extracted_seq_dict.values()), self._fasta_extracted)

        return concatenated_seqs
