# Functions to parse some of the metrics files from the WGS pipeline to generate the QC summary file

def parse_metric_file(metric_file, header_str, metric):
    """
    INPUT:
    - metric_file: path to the metric file
    - header_str:
        Picard metrics files contain the following information, preceeded by a header line:
        - Command
        - Date/time stamp
        - Summary Table containing the metrics
        - Histogram
        We need information from the summary table. To get the index of the first line of this table, we look forthe index
        of the specific header line (header_str).
    - metric: Exact metric we want to read from the summary table

    OUTPUT:
    Returns the metric
    """
    with open(metric_file) as handle:
        lines = handle.readlines()
        metrics_class_index = lines.index(header_str)

    header = lines[metrics_class_index + 1].split("\t")
    values = lines[metrics_class_index + 2].split("\t")

    parsed_metrics = dict(zip(header, values))

    return parsed_metrics[metric]

def get_fold_80(coverage_file, mean_coverage, percentile_20):
    """
    Fold80 is a way to calculate the evenness of coverage. It expresses the amount of additional sequencing needed to have
    80% of all targets covered at the currently observed mean. It is computed as the mean coverage divided by the 20th
    percentile.

    INPUT
    - coverage file: Output of Picard Raw WGS metrics file
    - Mean coverage: Mean coverage, parsed from the Raw WGS metrics file
    - Percentile_20: 20% of the total bases (genome_territory)

    COVERAGE FILE
      ## htsjdk.samtools.metrics.StringHeader
      # CollectRawWgsMetrics INPUT=/scratch/humangenomics/UZLeuven_run_WGS/GC088815/alignment/gather_bqsr_sorted_bam/gathered_sorted.bam OUTPUT=GC088815.raw.wgs.metrics.txt INCLUDE_BQ_HISTOGRAM=true INTERVALS=/db/refgenomes/Human/GATK-BroadIns/hg38/v0/wgs_coverage_regions.hg38.interval_list USE_FAST_ALGORITHM=true READ_LENGTH=250 TMP_DIR=[/temp/picard] VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=/db/refgenomes/Human/GATK-BroadIns/hg38/v0/Homo_sapiens_assembly38.fasta    MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=3 COVERAGE_CAP=100000 LOCUS_ACCUMULATION_CAP=200000 STOP_AFTER=-1 COUNT_UNPAIRED=false SAMPLE_SIZE=10000 ALLELE_FRACTION=[0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5] VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
      ## htsjdk.samtools.metrics.StringHeader
      # Started on: Wed Oct 06 12:11:19 CEST 2021

      ## METRICS CLASS        picard.analysis.CollectRawWgsMetrics$RawWgsMetrics
      GENOME_TERRITORY        MEAN_COVERAGE   SD_COVERAGE     MEDIAN_COVERAGE MAD_COVERAGE    PCT_EXC_ADAPTER PCT_EXC_MAPQ    PCT_EXC_DUPE    PCT_EXC_UNPAIRED        PCT_EXC_BASEQ   PCT_EXC_OVERLAP PCT_EXC_CAPPED  PCT_EXC_TOTAL    PCT_1X  PCT_5X  PCT_10X PCT_15X PCT_20X PCT_25X PCT_30X PCT_40X PCT_50X PCT_60X PCT_70X PCT_80X PCT_90X PCT_100X        HET_SNP_SENSITIVITY     HET_SNP_Q
      2745186691      49.838404       131.234792      48      5       0.000009        0       0.028816        0.000439        0.000004        0.08963 0.00041 0.119308        0.998591        0.996662        0.99479 0.9928980.99047  0.986888        0.980262        0.876739        0.439708        0.083201        0.014662        0.008625        0.006698        0.005369        0.996945        25

      ## HISTOGRAM    java.lang.Integer
      coverage        high_quality_coverage_count     unfiltered_baseq_count
      0       3868457 0
      1       1674733 0
      2       1310667 0
      [...]
      99999   0       0
      100000  881     0

    """
    with open(coverage_file) as handle:
        lines = handle.readlines()
        histogram_start_index = lines.index("## HISTOGRAM\tjava.lang.Integer\n")

    # Read the entire coverage histogram.
    ## columns: coverage - coverage_count - unfiltered baseq count
    coverage_histo = lines[histogram_start_index+2:-1]
    coverage_histo = [l.rstrip().split("\t") for l in coverage_histo]

    # get the list of coverage values. As we do not change the order, the index is the coverage
    # e.g. [7, 100, 500, 100] -> 7 bases 0X, 100 bases 1X, 500 bases 2X, 100 bases 3X
    cov_values = [int(l[1]) for l in coverage_histo]

    # calculate cumsum for the histogram:
    # e.g. [7, 100, 500, 100] -> [7, 107, 607, 707]: 7 bases 0X, 107 bases max. 1X, 607 bases max. 2X, 707 bases max. 3X
    histo_cumsum = [sum(cov_values[:i+1]) for i, value in enumerate(cov_values)]

    # add the input percentile_20 to the cumsum list and sort
    tmp = histo_cumsum + [int(percentile_20)]
    tmp.sort()
    # Find the coverage of the 20th percentile
    # e.g. 20% 707 = 141 -> [7, 107, 141, 607, 707]
    find_index = tmp.index(int(percentile_20))
    val_before = histo_cumsum[find_index - 1]
    val_after = histo_cumsum[find_index]
    # in our example, find_index will be 2, val_before: 107, val_after: 607

    ## Calculate the coverage for the lowest 20% of bases
    # val_before: 107: 107/707 = 15.13% <= 1X coverage
    # val_after: 607: 607/707 = 85.85% <= 2X coverage
    # 20% is between 15-85%: coverage of 20th percentile is between 1X and 2X
    # Fraction: (141-107) / (607-107+1) = 0.068
    # This is assuming linearity on the small interval, which is expected if the coverage follows a continuous distribution
    cov_20 = (((percentile_20 - val_before) / (val_after - val_before + 1)) + find_index - 1)

    # ideally, as close as possible to 1
    return mean_coverage / cov_20
