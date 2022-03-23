<!-- PROJECT TITLE -->
# Snakemake for Calling Variants

<!-- TABLE OF CONTENTS -->
## Table of contents
* [Introduction](#introduction)
* [Quick start](#quick-start)
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
files. Make sure you have enough storage space to run these scripts. These workflows are also suited for 
low-throughput DNA input data. DeepVariant workflow uses gpu, it is recommended to make sure gpu drivers are installed 
to avoid errors (see '[Running on a machine with GPU](https://github.com/google/deepvariant/blob/r0.9/docs/deepvariant-case-study.md)')
for details. 

<!-- Quick start -->
## Quick start 

<!-- Funding -->
## Funding
