from pathlib import Path

from camel.app.core.snakemake import snakemakeutils
from camel.scripts.broadwgs.snakefile import qc_summary


rule generate_qc_summary:
    """
    Generate a summary table containing all WGS QC metrics selected by WG4 from the 1+MG project (16/06/2021).
    """
    input:
        mark_duplicates_metrics = Path(config['working_dir']) / 'qc' / 'mark_duplicates' / "duplicate_metrics.txt.io",
        alignment_summary_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.alignment_summary_metrics',
        insert_size_metrics = Path(config['working_dir']) / "qc" / "aggregation_metrics" / f'{config["sample"]}.agg.insert_size_metrics',
        wgs_metrics = Path(config['working_dir']) / "qc" / "wgs_metrics" / f'{config["sample"]}.wgs.metrics.txt',
        variant_calling_metrics = Path(config['working_dir']) / "qc" / "variant_calling_metrics" / f'{config["sample"]}.variant_calling_metrics.io'
    output:
        QC_summary = Path(config['working_dir']) / "qc" / "QC_summary.txt"
    run:
        # Initialize dictionary containing threshold information as decided in EU 1+MG WG4
        alignment_metrics_threshold = {
            'percent_duplication' : {'val' : 0.1, 'threshold' : 'upper'},
            'percent_chimeras' : {'val': 0.01, 'threshold' : 'upper'},
            'median_insert_size' : {'val': 300, 'threshold': 'lower'},
            'mean_coverage' : {'val': 30, 'threshold': 'lower'},
            'evenness_coverage' : {'val': 1, 'threshold' : 'both'}
        }

        variant_calling_metrics_threshold = {
            "TITV_DBSNP" : {'val': 2, 'threshold': 'both'},
            "TITV_NOVEL" : {'val': 2, 'threshold': 'both'}
        }

        # initialize dictionaries to collect metrics
        alignment_metrics = {}
        variant_calling_metrics = {}

        ## Parse files and get relevant metrics
        # Mark Duplicates metrics
        mark_duplicates_file = Path(snakemakeutils.load_object(Path(input.mark_duplicates_metrics))[0].path)
        alignment_metrics['percent_duplication'] = qc_summary.parse_metric_file(mark_duplicates_file,
                                                                                "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n",
                                                                                "PERCENT_DUPLICATION")

        # Alignment summary metrics
        alignment_metrics['percent_chimeras'] = qc_summary.parse_metric_file(input.alignment_summary_metrics,
                                                                             "## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\n",
                                                                             'PCT_CHIMERAS')

        # Insert size metrics
        alignment_metrics['median_insert_size'] = qc_summary.parse_metric_file(input.insert_size_metrics,
                                                                               "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n",
                                                                               "MEDIAN_INSERT_SIZE")

        # WGS metrics
        mean_coverage = float(qc_summary.parse_metric_file(input.wgs_metrics,
                                                           "## METRICS CLASS\tpicard.analysis.WgsMetrics\n",
                                                           'MEAN_COVERAGE'))
        # total genome size
        genome_territory = int(qc_summary.parse_metric_file(input.wgs_metrics,
                                                           "## METRICS CLASS\tpicard.analysis.WgsMetrics\n",
                                                           'GENOME_TERRITORY'))
        # 20% of the genome size
        percentile_20 = genome_territory * 0.2

        alignment_metrics['mean_coverage'] = mean_coverage
        alignment_metrics['evenness_coverage'] = qc_summary.get_fold_80(input.wgs_metrics, mean_coverage, percentile_20)

        # Variant calling metrics
        var_calling_file = Path(snakemakeutils.load_object(Path(input.variant_calling_metrics))[0].path)
        variant_calling_metrics["TITV_DBSNP"] = qc_summary.parse_metric_file(var_calling_file,
                                                                             "## METRICS CLASS\tpicard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics\n",
                                                                             'DBSNP_TITV')
        variant_calling_metrics["TITV_NOVEL"] = qc_summary.parse_metric_file(var_calling_file,
                                                                             "## METRICS CLASS\tpicard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics\n",
                                                                             'NOVEL_TITV')

        # Print all alignment and variant calling metrics to output file
        alignment_metrics_metrics_names = ['mean_coverage', 'percent_duplication', 'median_insert_size', 'percent_chimeras', 'evenness_coverage']
        variant_calling_metrics_names = ['TITV_DBSNP', 'TITV_NOVEL']
        with open(output.QC_summary, 'w') as outhandle:
            for aln_metric in alignment_metrics_metrics_names:
                outhandle.write("\t".join([aln_metric, str(alignment_metrics[aln_metric]),
                                           str(alignment_metrics_threshold[aln_metric]['val']),
                                           alignment_metrics_threshold[aln_metric]['threshold']]) + "\n")

            for var_metric in variant_calling_metrics_names:
                outhandle.write("\t".join([var_metric, str(variant_calling_metrics[var_metric]),
                                           str(variant_calling_metrics_threshold[var_metric]['val']),
                                           variant_calling_metrics_threshold[var_metric]['threshold']]) + "\n")

