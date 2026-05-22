import re
from typing import Optional

from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.picard.picard import Picard


class CollectMultipleMetrics(Picard):
    """
    Class for Picard CollectMultipleMetrics function to calculate various QC metrics for read mapping

    PROGRAM: Set of metrics programs to apply during the pass through the SAM file.
        Default values:
        - CollectAlignmentSummaryMetrics
        - CollectBaseDistributionByCycle
        - CollectInsertSizeMetrics
        - MeanQualityByCycle
        - QualityScoreDistribution

        Possible values:
        - CollectAlignmentSummaryMetrics
        - CollectInsertSizeMetrics
        - QualityScoreDistribution
        - MeanQualityByCycle
        - CollectBaseDistributionByCycle
        - CollectGcBiasMetrics
        - RnaSeqMetrics
        - CollectSequencingArtifactMetrics
        - CollectQualityYieldMetrics

    LEVEL/METRIC_ACCUMULATION_LEVEL: The level(s) at which to accumulate metrics
        Default value:
        - ALL_READS

        Possible values:
        - ALL_READS
        - SAMPLE
        - LIBRARY
        - READ_GROUP

    The PROGRAM and LEVEL options can be set to 'null' to clear the default list
    Usage example:
        java -jar picard.jar CollectMultipleMetrics \
        I=input.bam \
        O=multiple_metrics \
        R=reference_sequence.fasta \
        PROGRAM=null \
        PROGRAM=QualityScoreDistribution \
        PROGRAM=MeanQualityByCycle \
        LEVEL=null \
        LEVEL=LIBRARY
    """
    # Suffices for all PROGRAM options
    OUTPUT_FILE_SUFFIX = {
        'CollectAlignmentSummaryMetrics': {'AlignmentSummary': '.alignment_summary_metrics'},
        'CollectInsertSizeMetrics': {'InsertSize': '.insert_size_metrics', 'InsertSizeFigure': '.insert_size_histogram.pdf'},
        'QualityScoreDistribution': {'QualityDistribution': '.quality_distribution_metrics', 'QualityDistributionFigure': '.quality_distribution.pdf'},
        'MeanQualityByCycle': {'QualityByCycle': '.quality_by_cycle_metrics', 'QualityByCycleFigure': '.quality_by_cycle.pdf'},
        'CollectBaseDistributionByCycle': {'BaseDistributionByCycle': '.base_distribution_by_cycle_metrics', 'BaseDistributionByCycleFigure': '.base_distribution_by_cycle.pdf'},
        'CollectGcBiasMetrics': {'GcBias': '.gc_bias.detail_metrics', 'GcBiasSummary': '.gc_bias.summary_metrics', 'GcBiasFigure': '.gc_bias.pdf'},
        'RnaSeqMetrics': {'RnaSeq': '.rna_metrics'},
        'CollectSequencingArtifactMetrics': {'SequencingArtefactDetail': '.bait_bias_detail_metrics', 'SequencingArtefactSummary': '.bait_bias_summary_metrics',
                                             'SequencingArtefactErrorSummary': '.error_summary_metrics', 'SequencingArtefactPreAdapterDetail': '.pre_adapter_detail_metrics',
                                             'SequencingArtefactPreAdapterSummary': '.pre_adapter_summary_metrics'},
        'CollectQualityYieldMetrics': {'QualityYield': '.quality_yield_metrics'}
    }

    def __init__(self):
        """
        Initialize a picard tool
        :return: None
        """
        super().__init__('Picard CollectMultipleMetrics', '2.23.3')

        self._required_inputs = ['BAM', 'SAM']
        self._specific_parameters = ['metric_accumulation_level_multi']

    def _set_output(self) -> None:
        """
        Set the output for Picard CollectMultipleMetrics function
        :return: None
        """
        for key in self._parameters:
            if re.match('metrics_', key):
                # strip 'metrics_' from key
                tool = key.split("_")[1]

                # loop over all possible output suffices for that specific tool
                try:
                    for suffix_key, suffix_value in CollectMultipleMetrics.OUTPUT_FILE_SUFFIX[tool].items():
                        output_file = self.folder / (self._parameters['output_prefix'].value + suffix_value)
                        self._tool_outputs['TXT_' + suffix_key] = [ToolIOFile(output_file)]
                except KeyError:
                    logger.warning(f'Picard CollectMultipleMetrics unsupported metrics {key}, its results will not be analyzed or returned.')

    def _set_informs(self, stderr: str | None = None) -> None:
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        # Known Metrics: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,
        #   MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics,
        #   CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
        for key in self._parameters:
            if re.match('metrics_', key):
                if key == 'metrics_CollectAlignmentSummaryMetrics':
                    self.__analyze_alignment_summary()
                elif key == 'metrics_CollectInsertSizeMetrics':
                    self.__analyze_insert_size()
                elif key == 'metrics_CollectGcBiasMetrics':
                    self.__analyze_gc_bias()
                elif key == 'metrics_QualityScoreDistribution':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                elif key == 'metrics_CollectBaseDistributionByCycle':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                elif key == 'metrics_MeanQualityByCycle':
                    # No output to analyze, set to prevent going to 'else' case
                    continue
                else:  # RnaSeqMetrics, CollectSequencingArtifactMetrics, CollectQualityYieldMetrics
                    logger.warning(
                        f'Picard CollectMultipleMetrics unsupported metrics {key}, its results will not be analyzed or returned.')

    def _build_command(self, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Build the command to run tool
        :return: None
        """
        build_options = self._build_options(excluded_parameters=self._specific_parameters, delimiter='=')

        if 'metric_accumulation_level_multi' in self._parameters:
            build_options.append(self.__split_multi_options('metric_accumulation_level_multi'))

        option_string = " ".join(build_options)

        self._command.command = " ".join([
            *self._get_base_command(), self._java_options_temp_dir,
            self._input_string, self._output_string, option_string
        ])

    def __analyze_alignment_summary(self) -> None:
        """
        Analyze the Alignment Summary statistics
        :return: None
        """
        # Ref: http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics
        align_stats = {}
        with open(self._tool_outputs['TXT_AlignmentSummary'][0].path) as inf:
            col_nbs = None
            for line in inf:
                if not line.startswith('#') and 'TOTAL_READS' in line and col_nbs is None:  # First line after comments is the header
                    col_nbs = {x: i for i, x in enumerate(line.strip().split('\t'))}
                if re.match('PAIR', line):
                    pair_informs = line.split('\t')
                    align_stats['TOTAL_READS'] = pair_informs[col_nbs['TOTAL_READS']]
                    # PF: Passing (Illumina) Filter
                    align_stats['PassingFilter_READS'] = pair_informs[col_nbs['PF_READS']]
                    align_stats['PF_NOISE_READS'] = pair_informs[col_nbs['PF_NOISE_READS']]
                    align_stats['PF_READS_ALIGNED'] = pair_informs[col_nbs['PF_READS_ALIGNED']]
                    align_stats['PCT_PF_READS_ALIGNED'] = pair_informs[col_nbs['PCT_PF_READS_ALIGNED']]
                    align_stats['PF_ALIGNED_BASES'] = pair_informs[col_nbs['PF_ALIGNED_BASES']]
                    align_stats['PF_MISMATCH_RATE'] = pair_informs[col_nbs['PF_MISMATCH_RATE']]
                    align_stats['PF_INDEL_RATE'] = pair_informs[col_nbs['PF_INDEL_RATE']]
                    align_stats['MEAN_READ_LENGTH'] = pair_informs[col_nbs['MEAN_READ_LENGTH']]
                    align_stats['READS_ALIGNED_IN_PAIRS'] = pair_informs[col_nbs['READS_ALIGNED_IN_PAIRS']]
                    align_stats['PCT_READS_ALIGNED_IN_PAIRS'] = pair_informs[col_nbs['PCT_READS_ALIGNED_IN_PAIRS']]

                    break

        self.informs['AlignmentSummary_stats'] = align_stats

    def __analyze_insert_size(self) -> None:
        """
        Analyze the Insert Size statistics
        :return: None
        """
        # Ref: http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics
        data_row = False
        self.informs['InsertSize_stats'] = {}
        with open(self._tool_outputs['TXT_InsertSize'][0].path) as inf:
            for line in inf:
                if re.match('MEDIAN_INSERT_SIZE', line):
                    data_row = True
                    col_nbs = {x: i for i, x in enumerate(line.strip().split('\t'))}
                    continue
                if data_row:
                    if line.strip() == '':
                        break
                    informs = line.split('\t')
                    # Only take statistics for FR (READ_PAIRS) reads
                    if informs[col_nbs['PAIR_ORIENTATION']] == 'FR':
                        self.informs['InsertSize_stats']['MEDIAN_INSERT_SIZE'] = informs[col_nbs['MEDIAN_INSERT_SIZE']]
                        self.informs['InsertSize_stats']['MEDIAN_ABSOLUTE_DEVIATION'] = informs[col_nbs['MEDIAN_ABSOLUTE_DEVIATION']]
                        self.informs['InsertSize_stats']['MEAN_INSERT_SIZE'] = informs[col_nbs['MEAN_INSERT_SIZE']]
                        self.informs['InsertSize_stats']['STANDARD_DEVIATION'] = informs[col_nbs['STANDARD_DEVIATION']]
                    else:
                        logger.warning(
                            f'Picard CollectMultipleMetrics unsupported READS ORIENTATION {informs[col_nbs["PAIR_ORIENTATION"]]} for InsertSize analysis')

    def __analyze_gc_bias(self) -> None:
        """
        Analyze the GC Bias statistics
        :return: None
        """
        # Only GC Bias summary is analyzed
        self.informs['GCBias_stats'] = {}
        with open(self._tool_outputs['TXT_GcBiasSummary'][0].path) as inf:
            for line in inf:
                if re.match('All Reads', line):
                    informs = line.split("\t")
                    self.informs['GCBias_stats']['AT_DROPOUT'] = informs[4]
                    self.informs['GCBias_stats']['GC_DROPOUT'] = informs[5]
                    break

    def __split_multi_options(self, option) -> str:
        """
        Multiple values allowed for certain parameters. These are passed in a comma separated string and need to be split
        """
        option_list = self._parameters[option].value.split(",")
        return "".join(f" {self._parameters[option].option}={s} " for s in option_list)
