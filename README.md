<!-- PROJECT TITLE -->
# Snakemake for Calling Variants

<!-- TABLE OF CONTENTS -->
## Table of contents
* [Introduction](#introduction)
* [Quick start](#quick-start)
* [Debugging](#debugging)
* [Paper and Citation](#paper-and-citation)
* [Funding](#funding)

<!-- Introduction -->
## Introduction
The snakemake scripts are used to automate call variants given one or more pairs of fastq files. There are two
approaches available to call variants. 
1. Use bcftools
2. Use DeepVariant

Both approaches use a configuration file that makes the workflow more flexible and abstracts away direct dependencies. 
Data pre-processing steps such as indel realignment, mark duplicates, and BQSR are part of the workflow. However, DeepVariant 
workflow doesn't include BQSR as it decreases the accuracy of variant calling. These workflows create temporary intermediate 
files. Make sure to have enough storage space to run these scripts. These workflows are also suited for 
low-throughput DNA input data. DeepVariant workflow uses gpu, it is recommended to make sure gpu drivers are installed 
to avoid errors (see '[Running on a machine with GPU](https://github.com/google/deepvariant/blob/r0.9/docs/deepvariant-case-study.md)')
for details. For quick snakemake example see [this](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html).

<!-- Quick start -->
## Quick start
* This assumes latest version of snakemake, bwa-mem, bcftools, deepvariant, and gatk tools are installed. 

* Provide the directory where the reads are located. Provide only the read 1 fastq file as input. 
If you have multiple samples, provide a sample specific regex pattern for read 1. Illumina, by default, uses underscores in
their naming convention. The snakemake scripts work on these underscores to extract sample name and get its corresponding 
read 2. The sample name should not have any underscore in it, dashes and dots should work fine.

* Change config file to include the correct paths to grch38, dbsnp_138.hg38.vcf.gz, Mills_and_1000G_gold_standard.indels.hg38.vcf.gz, 
and 1000G_phase1.snps.high_confidence.hg38.vcf.gz and their indices.

* In snakemake, point to the correct picard.jar.

* Make sure gatk is added to the environment variable.

* Make sure correct path to tmpdir is mentioned in snakemake, markdup_bams. 

* Bcftools snakemake: 
```bash
snakemake -s snakemake_picard_bcftools.smk --config Reads="example/Sample1*R1*.fastq.gz" -c128
```

* Bcftools deepvariant:
```bash
snakemake -s snakemake_picard_deepvariant.smk --config Reads="example/Sample1*R1*.fastq.gz" -c128
```

* To just see the steps that will be run, it is recommended to do a dry run. For example,
```shell
snakemake -nps snakemake_picard_bcftools.smk --config Reads="example/Sample1*R1*.fastq.gz" -c128
```

<!--Debugging-->
## Debugging
If you get errors or snakemake stops after any particular step, try running the same command by copying and
pasting it. Sometimes, just running the command by itself gives you more details about the error. 

<!--Paper and Citation-->
## Paper and Citation
Woerner AE, Mandape S, Kapema KB, Duque TM, Smuts A, King JL, Crysup B, Wang X, Huang M, Ge J, Budowle B. Optimized Variant Calling 
for Estimating Kinship. Forensic Science International: Genetics (2022); 61(102785)

<!--Funding-->
## Funding
This work was supported in part by award 2019-DU-BX-0046, awarded by the National Institute of Justice, Office of
Justice Programs, U.S. Department of Justice and by internal funds from the Center for Human
Identification. The opinions, findings, and conclusions or recommendations expressed in this publication are those of 
the authors and do not necessarily reflect those of the U.S. Department of Justice.
