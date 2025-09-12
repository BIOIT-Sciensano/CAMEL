# Overview
The *Neisseria* pipeline performs complete characterization of *Neisseria meningitidis* isolates.

Version: **1.4**

# Components

**Note:** If the input type is `fasta`, pre-processing steps 1 to 4 are skipped.

## 1. Human read removal (optional) 

If enabled, human reads are removed using the NCBI Human Read Removal Tool (HRRT) 2.2.1.
The tool is executed with default options.

## 2. Coverage check

The workflow starts by checking the coverage of the input FASTQ datasets. 
Coverage is estimated by dividing the total number of bases by the size of the `NZ_CP021520.1` *N. meningitidis* 
11-7 reference genome. The total number of bases in the FASTQ files is determined using the `size` function of 
`seqtk 1.4`.

Datasets with an estimated coverage >=100x are downsampled to ~100x using the `subsample` function of `seqtk 1.4`.

## 3. Read trimming

### Illumina

Read trimming is performed using `fastp 0.23.4` (default) or `trimmomatic 0.39`.

For `fastp` the following options are used:
```
--compression 4
--detect_adapter_for_pe
--cut_front
--cut_front_window_size 1
--cut_front_mean_quality 10
--cut_tail
--cut_tail_window_size 1
--cut_tail_mean_quality 10
--cut_right
--cut_right_window_size 4
--cut_right_mean_quality 20
--length_required 40
```

For `trimmomatic` the following options are used:
```
-phred33
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10
LEADING:10 
TRAILING:10 
SLIDINGWINDOW:4:20 
MINLEN:40
```

Quality reports are generated before and after trimming using `fastqc 0.12.1`.

### ONT

Read filtering is performed using `seqkit 2.3.1` using the following options:

```
--min-len 1000
--min-qual 10
```
**Note:** these values can be changed using the command-line options.

Quality reports are generated before and after filtering using `NanoPlot 1.41.6`.

## 4. Assembly

### Illumina

Processed reads are assembled using `SPAdes 3.15.5` with the following options:
```
--cov-cutoff 10
--isolate
```

### ONT

Filtered reads are assembled using `Flye 2.9.4` with default options providing the filtered reads using the 
`--nano-corr` option.

### QUAST

`QUAST 5.2.0` is then used to check the quality of the resulting assembly with the following options:

```
-r {ref_genome_fasta}
--features {ref_genome_gff3}
--pe1 {forward_reads}
--pe2 {reverse_reads}
```

For ONT input, the `--pe1` and `--pe2` options are replaced by `--nanopore`.

For FASTA input, these options are omitted.

### BUSCO

The completeness of the assembly is checked using `BUSCO 5.5.0` with the following options:
```
--mode genome
--offline
--lineage_dataset bacteria_odb10
```

## 5. Advanced QC

### Kraken 2

The trimmed reads or contigs are checked for contamination using `kraken2 2.1.1` against an in-house database with 
microbial genomes. 
The date of the last database update is included in the output report.

### ConFindr

The samples are screened for inter- and intra-species contamination using `ConFindr 0.8.2` with the ribosomal MLST 
database.

**Note:** ConFindr is not executed when the input type is `fasta`.

**Note:** ConFindr on ONT input data is experimental.

### Quality checks

An overview of the quality checks is provided below. Warnings are included for quality checks that fail but do not stop 
the pipeline execution. 


| **metric**                             | **warning threshold** | **fail threshold** | **description**                                                                                                                                                                                                                 |
|----------------------------------------|-----------------------|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kraken: contaminants                   | 1.00%                 | 5.00%              | Percentage of reads / contigs assigned to species other than *N. meningitidis*                                                                                                                                                  |
| Typing loci detected (%)               | 90%                   | 95%                | Percentage of cgMLST loci detected (or MLST loci when cgMLST is disabled)                                                                                                                                                       |
| Coverage against assembled contigs     | 20x                   | 10x                | Coverage of the reads mapped to the assembly (determined by QUAST)                                                                                                                                                              |
| Reads mapping to the assembled contigs | 95%                   | 90%                | Percentage of reads mapping back to the assembly (determined by QUAST)                                                                                                                                                          |
| Total assembly length deviation        | 10%                   | 20%                | Percent deviation from the expected genome size (determined from the reference genome)                                                                                                                                          |
| ConFindr: number of contaminating SNPs | 10                    | 20                 | Number of SNPs flagged as contaminant by ConFindr                                                                                                                                                                               | 
| Percentage of complete BUSCO genes     | 90%                   | 95%                | Percentage of complete BUSCO genes identified                                                                                                                                                                                   |
| FastQC: Average quality score          | 30                    | 25                 | Checks if the average read quality is above the given threshold.                                                                                                                                                                |
| FastQC: GC-content deviation           | 2.00%                 | 4.00%              | Checks if the detected GC content is close enough to the expected GC content for this organism (38.00%).                                                                                                                        |
| FastQC: Max. N-fraction                | 0.0050                | 0.0100             | Checks if the maximal N fraction at any read position is below the given threshold.                                                                                                                                             |
| FastQC: Per-base sequence content      | 3.00%                 | 6.00%              | Checks if the difference between A-T and C-G is below the given threshold at every position. The first 20 and last 5 bases of the reads are skipped, as the peaks there can be caused by the library kit or trimming artifacts. |
| FastQC: Q-score drop                   | 200                   | 150                | Checks whether the average position in the reads where the mean Q-score drops below 30 is above the given threshold.                                                                                                            |
| FastQC: Sequence length distribution   | 66.67%                | 40.00%             | Checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |
| NanoPlot: Median read length           | 500                   | 250                | Checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |
| NanoPlot: Median read quality          | 10                    | 8                  | Checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |
| seqkit: GC-content deviation           | 2.00%                 | 4.00%              | Checks if the median read length of the trimmed reads is below a threshold compared to the mode length of the raw input reads (251).                                                                                            |

**Note:** FastQC metrics are evaluated separately for the forward and reverse reads.

The QC checks enabled for the supported input types are listed in the table below.

| **metric**                             | **illumina** | **fasta** | **ont** |
|----------------------------------------|--------------|-----------|---------|
| Kraken: contaminants                   | Yes          | Yes       | Yes     |
| Typing loci detected (%)               | Yes          | Yes       | Yes     |  
| Coverage against assembled contigs     | Yes          | No        | Yes     |
| Reads mapping to the assembled contigs | Yes          | No        | Yes     |
| Total assembly length deviation        | Yes          | Yes       | Yes     |
| ConFindr: number of contaminating SNPs | Yes          | No        | Yes     | 
| Percentage of complete BUSCO genes     | Yes          | Yes       | Yes     |
| FastQC: Average quality score          | Yes          | No        | No      |
| FastQC: GC-content deviation           | Yes          | No        | No      |
| FastQC: Max. N-fraction                | Yes          | No        | No      |
| FastQC: Per-base sequence content      | Yes          | No        | No      |
| FastQC: Q-score drop                   | Yes          | No        | No      |
| FastQC: Sequence length distribution   | Yes          | No        | No      |
| NanoPlot: Median read length           | No           | No        | Yes     |
| NanoPlot: Median read quality          | No           | No        | Yes     |
| seqkit: GC-content deviation           | No           | No        | Yes     |

## 6. Variant calling & filtering (optional)

Reads are mapped against the reference genome using `Bowtie2 2.5.1`, for illumina data and `minimap2 2.26` for ONT data.
Variants are then called with `bcftools 1.17`, using the `bcftools mpileup` followed by `bcftools call` with the 
following options:

```
bcftools mpileup {BAM in} --fasta-ref {FASTA ref} --output-type z --count-orphans --output {PILEUP out};
bcftools call {PILEUP out} --consensus-caller --output {VCF_GZ out} --output-type z --variants-only --ploidy 1;
```

The following variant filters then applied, with threshold values listed in the output report.

- Depth (see command in the output report)
- SNP/indel quality (see command in the output report)
- Mapping quality (see command in the output report)
- Distance (in-house script to filter SNPs located within 10 bp of another SNP)
- Z-score (in-house script to filter based on Z-score & Y-multiplier as described by [Kaas et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4128722/))

## 7. Antimicrobial resistance (AMR) characterization

### NCBI AMRFinder+

`NCBI AMRFinder+ 3.11.26` is used with default options to detect genes and mutation associated with AMR.

The `--organism` option is set to 'Neisseria meningitidis'

The database version is indicated in the output report and summary output file.

### ResFinder4

`ResFinder4 4.4.2` is used with the following options to detect genes and mutation associated with AMR:

```
--min_cov 0.6
--acquired
--threshold 0.9
```

The database version is indicated in the output report and summary output file.

## 8. Sequence typing

Sequence typing is performed as described in [Bogaerts *et al.*](https://pubmed.ncbi.nlm.nih.gov/30894839/) with an 
updated version of blast (`blast 2.14.0`). 
Alternative detection using `kma 1.4.12a` or `SRST2 0.2.0` is available by changing the `--detection-method` parameter.

The following typing schemes are available:

| **name**                 | **origin** |
|--------------------------|------------|
| rMLST                    | PubMLST    |
| Classic MLST             | PubMLST    |
| rplF                     | PubMLST    |
| PorA                     | PubMLST    |
| PorB                     | PubMLST    |
| FetA                     | PubMLST    |
| Antibiotic resistance    | PubMLST    |
| Vaccine targets          | PubMLST    |
| Factor-H binding protein | PubMLST    |
| cgMLST                   | PubMLST    |

## 9. Antigen typing

### Bexsero antigen sequence typing (BAST)

Antigen typing is based on typing using the BAST scheme from PubMLST, using the method described in the 
*Sequence typing* section.

### gMATS

gMATS is used to predict the efficacy of the Bexsero vaccine. It works by matching the alleles of the BAST typing 
scheme to a database collected from [literature](https://doi.org/10.1016/j.vaccine.2018.12.061).

## 10. Serotype determination

Serotype determination is performed using the [characterize_neisseria_capsule](https://github.com/ntopaz/characterize_neisseria_capsule) script (commit `a75a009`).
