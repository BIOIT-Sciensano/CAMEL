import dataclasses
import operator
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Any

SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_REPORT = 'quality_checks/report/html.iob'
OUTPUT_SUMMARY = 'quality_checks/summary/summary_out.{ext}'
OUTPUT_REPORT_JSON = 'quality_checks/report/informs.json'


@dataclasses.dataclass
class QCCheck:
    """
    Dataclass to represent quality checks.
    """
    key: str
    full_name: str
    threshold_warn: float
    threshold_fail: float
    supported_input_types: list = dataclasses.field(default_factory=list)
    fmt_string_value: Optional[str] = None
    value_should_exceed: Optional[bool] = True
    explanation: Optional[str] = None
    is_default: bool = True

    def to_dict(self, value: Optional[float] = None) -> dict[Any, str]:
        """
        Converts the QC check to a dictionary.
        :param value: Metric value
        :return: Dictionary
        """
        data = asdict(self)
        comp = operator.lt if self.value_should_exceed else operator.gt
        if value is None:
            data['status'] = 'Skipped'
        elif comp(value, self.threshold_fail):
            data['status'] = 'Failed'
        elif comp(value, self.threshold_warn):
            data['status'] = 'Warning'
        else:
            data['status'] = 'OK'
        data['value'] = value
        return data


QC_CHECKS_BY_KEY = {qc.key: qc for qc in [
    QCCheck(
        key='kraken_illumina',
        full_name='Kraken 2: contaminants (Illumina)',
        threshold_warn=1.0,
        threshold_fail=5.0,
        fmt_string_value='{:.2f}%',
        supported_input_types=['hybrid', 'illumina'],
        value_should_exceed=False),
    QCCheck(
        key='kraken_ont',
        full_name='Kraken 2: contaminants (ONT)',
        threshold_warn=1.0,
        threshold_fail=5.0,
        fmt_string_value='{:.2f}%',
        supported_input_types=['hybrid', 'ont'],
        value_should_exceed=False),
    QCCheck(
        key='kraken_fasta',
        full_name='Kraken 2: contaminants (FASTA)',
        threshold_warn=1.0,
        threshold_fail=5.0,
        fmt_string_value='{:.2f}%',
        supported_input_types=['fasta'],
        value_should_exceed=False),
    QCCheck(
        key='confindr',
        full_name='ConFindr: number of contaminating SNPs',
        threshold_warn=25,
        threshold_fail=50,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:,}',
        value_should_exceed=False),
    QCCheck(
        key='typing_loci',
        full_name='Typing loci detected (%)',
        threshold_warn=95.0,
        threshold_fail=90.0,
        fmt_string_value='{:.2f}%',
        supported_input_types=['hybrid', 'illumina', 'ont', 'fasta', 'fasta_with_vcf']),
    QCCheck(
        key='map_rate_ref_illumina',
        full_name='Reads mapping to reference genome (Illumina)',
        threshold_warn=95.0,
        threshold_fail=90.0,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}%',
        is_default=False),
    QCCheck(
        key='map_rate_ref_ont',
        full_name='Reads mapping to reference genome (ONT)',
        threshold_warn=95.0,
        threshold_fail=90.0,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}%',
        is_default=False),
    QCCheck(
        key='map_rate_assembly_illumina',
        full_name='Reads mapping to the assembled contigs (Illumina)',
        threshold_warn=95.0,
        threshold_fail=90.0,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}%'),
    QCCheck(
        key='map_rate_assembly_ont',
        full_name='Reads mapping to the assembled contigs (ONT)',
        threshold_warn=95.0,
        threshold_fail=90.0,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}%'),
    QCCheck(
        key='cov_assembly_illumina',
        full_name='Coverage against the assembled contigs (Illumina)',
        threshold_warn=20,
        threshold_fail=10,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}x'),
    QCCheck(
        key='cov_ref_illumina',
        full_name='Coverage against the reference genome (Illumina)',
        threshold_warn=20,
        threshold_fail=10,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}x',
        is_default=False),
    QCCheck(
        key='cov_assembly_ont',
        full_name='Coverage against the assembled contigs (ONT)',
        threshold_warn=20,
        threshold_fail=10,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}x'),
    QCCheck(
        key='cov_ref_ont',
        full_name='Coverage against the reference genome (ONT)',
        threshold_warn=20,
        threshold_fail=10,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}x',
        is_default=False),
    QCCheck(
        key='assembly_total_len',
        full_name='Total assembly length deviation',
        threshold_warn=10,
        threshold_fail=20,
        supported_input_types=['fasta', 'fasta_with_vcf', 'hybrid', 'illumina', 'ont'],
        fmt_string_value='{:.2f}%',
        value_should_exceed=False,
        explanation='Percent deviation from the expected genome size.'),
    QCCheck(
        key='busco',
        full_name='Percentage of complete BUSCO genes',
        threshold_warn=90.0,
        threshold_fail=80.0,
        fmt_string_value='{:.2f}%',
        supported_input_types=['fasta', 'fasta_with_vcf', 'hybrid', 'illumina', 'ont'],
        explanation='BUSCO Benchmarking Universal Single-Copy Orthologs) can be used to estimate the completeness of '
                    'an assembly.'),
    QCCheck(
        key='fqc_avg_qual_{ori}',
        full_name='FastQC: Average quality score',
        threshold_warn=30.0,
        threshold_fail=25.0,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.0f}',
        explanation='checks if the average read quality is above the given threshold.'),
    QCCheck(
        key='fqc_gc_{ori}',
        full_name='FastQC: GC-content deviation',
        threshold_warn=2.0,
        threshold_fail=4.0,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}%',
        value_should_exceed=False,
        explanation='checks if the detected %GC-content is close enough to the expected value for this organism '
                    '(<b>{:.2f}%</b>).'),
    QCCheck(
        key='fqc_n_fraction_{ori}',
        full_name='FastQC: Max. N-fraction',
        threshold_warn=0.025,
        threshold_fail=0.050,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.4f}',
        value_should_exceed=False,
        explanation='checks if the maximal N fraction at any read position is below the given threshold.'),
    QCCheck(
        key='fqc_per_base_{ori}',
        full_name='FastQC: Per-base sequence content',
        threshold_warn=3.0,
        threshold_fail=6.0,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}%',
        value_should_exceed=False,
        explanation='checks if the difference between A-T and C-G is below the given threshold at every position. The '
                    'first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by the '
                    'library kit or trimming artifacts.'),
    QCCheck(
        key='fqc_qscore_{ori}',
        full_name='FastQC: Q-score drop',
        threshold_warn=200,
        threshold_fail=150,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.0f}',
        explanation='checks whether the average position in the reads where the mean Q-score drops below <b>30</b> is '
                    'above the given threshold.'),
    QCCheck(
        key='fqc_seq_len_{ori}',
        full_name='FastQC: Sequence length distribution',
        threshold_warn=66.67,
        threshold_fail=40.00,
        supported_input_types=['hybrid', 'illumina'],
        fmt_string_value='{:.2f}%',
        explanation='checks if the median read length of the trimmed reads is below a threshold compared to the mode '
                    'length of the raw input reads (<b>{}</b>).'),
    QCCheck(
        key='nanoplot_len',
        full_name='NanoPlot: Median read length (ONT)',
        threshold_warn=500,
        threshold_fail=250,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:,}'),
    QCCheck(
        key='nanoplot_qual',
        full_name='NanoPlot: Median read quality (ONT)',
        threshold_warn=10,
        threshold_fail=8,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}'),
    QCCheck(
        key='seqkit_gc',
        full_name='seqkit: GC-content deviation (ONT)',
        threshold_warn=2.0,
        threshold_fail=4.0,
        supported_input_types=['hybrid', 'ont'],
        fmt_string_value='{:.2f}%',
        value_should_exceed=False,
        explanation='checks if the detected GC content is close enough to the expected GC content for this organism'
                    ' (<b>{:.2f}%</b>).')
]}


def get_qc_checks(input_type: str, skipped_checks: list[str] = None, forced_checks: Optional[list[str]] = None) \
        -> list[Path]:
    """
    Returns the output paths for the QC checks of the corresponding input type.
    :param input_type: Input type
    :param skipped_checks: (Optional) list of skipped QC checks
    :param forced_checks: (Optional list of QC checks that are forced
    :return: List of paths with the enabled QC checks
    """
    paths = []
    if forced_checks is None:
        forced_checks = []

    # Add defaults
    for key, qc_check in QC_CHECKS_BY_KEY.items():
        if input_type not in qc_check.supported_input_types:
            continue
        elif key in skipped_checks:
            continue
        elif (qc_check.is_default is False) and (key not in forced_checks):
            continue

        # Add orientation for PE QC checks
        if '{ori}' in key:
            paths.extend([Path('quality_checks', f'{key.format(ori=ori)}.json') for ori in ('fwd', 'rev')])
        else:
            paths.append(Path('quality_checks', f'{key}.json'))

    # Add forced QC checks
    return paths
