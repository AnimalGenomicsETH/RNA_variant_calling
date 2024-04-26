[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![snakemaker](https://github.com/AnimalGenomicsETH/RNA_variant_calling/actions/workflows/snakemaker.yaml/badge.svg?branch=main)](https://github.com/AnimalGenomicsETH/RNA_variant_calling/actions/workflows/snakemaker.yaml)

## RNA sequencing variants are enriched for eQTL in cattle tissues

DNA sequencing is widely used for calling variants while RNA sequencing is relegated to only profiling gene expression.
Given recent advancements in [RNA sequencing variant calling](https://github.com/google/deepvariant/blob/r1.6.1/docs/deepvariant-rnaseq-case-study.md), we wanted to see how far we could push the RNA.

Broadly, total RNA sequencing in three different cattle tissues calls far more variants than we were expecting, with relatively good precision.
However, there are unresolved RNA DNA differences (RDDs) that mean we can't always trust the RNA-seq variants, even after **conservative** filtering.
This can be a problem when imputing RNA variants into large reference panels, as these RDDs or other RNA-specific effects like *allele specific expression* can interfer.
We even find some strange eQTL when using RNA-seq variants.

### Cite

Preprint is coming soon.

### Using the code

All the DNA- and RNA-seq data comes from [Mapel et al. 2024](https://www.nature.com/articles/s41467-024-44935-7), and these pipelines should be useable (up to maybe some specific features to our Euler cluster).
The main steps are as follows:
 - Variant calling from aligned BAMs
 - Assessing F1 score between DNA- and RNA-seq variants
 - Building gene expression matricies from RNA alignments
 - Association mapping between DNA- or RNA-seq variants and gene expression
