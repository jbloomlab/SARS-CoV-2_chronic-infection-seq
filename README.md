## SARS-CoV-2 Chronic Infection Sequencing
*2020-11-22*

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.17-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4433186.svg)](https://doi.org/10.5281/zenodo.4433186)

### Authors

* [Will Hannon](https://www.linkedin.com/in/williamhannon/)
* [Manish Choudhary](https://jonathanlilab.bwh.harvard.edu/)
* [Jonathan Li](https://jonathanlilab.bwh.harvard.edu/)
* [Jesse Bloom](https://www.fredhutch.org/en/faculty-lab-directory/bloom-jesse.html)

### Overview

This repo contains the code for generating the analysis of intra-host viral variation from the deep-sequencing data of SARS-CoV-2 from an immunocompromised patient. The associated publication for these samples is located [here](https://www.nejm.org/doi/full/10.1056/NEJMc2031364), and the preprint is located [here](https://www.biorxiv.org/content/10.1101/2020.11.30.405472v1).

Intra-patient single-nucleotide polymorphisms (SNPs) were identified with an automated variant-calling pipeline created with `Snakemake` (Köster & Rahmann, 2012). Briefly, paired-end reads were filtered, and sequencing adaptors were removed with `fastp` (Chen et al., 2018). Reads from SARS-CoV-2 were enriched by kmer matching to the Wuhan-Hu-1 reference genome (NC_045512.2) using [`BBDuk`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Following filtering, reads were aligned to the Wuhan-Hu-1 reference with `BWA-MEM` (Li, 2013). Variants were identified by counting the coverage of each base at every position in the reference genome using a custom Python script. These variants were filtered based on a minimum allele frequency of >0.01, a PHRED quality threshold of >25, and coverage of more than 100 reads.

 * Samples on the SRA have already been filtered using `fastp` to remove any human reads before uploading.
 
## Getting Started 

First, clone the repository to your desired location. 

```
git clone https://github.com/jbloomlab/SARS-CoV-2_chronic-infection-seq.git 
```

Once you've cloned the repository and installed conda, you can create the environment needed to run the pipeline. To do this, run the following command. 

```
conda env create --file environment.yml; conda activate viral-deepseq
```

## Running Analysis

To configure the analysis, you will need a table formatted like the one below. This table should contain the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) accessions, the Library Layout (paired-end or single-end files), if the files are single-ended format but interleaved, the name of the virus corresponding to the location of its genome in the [`config file`](/config/config.yml), and the identity of the host organism if there are contaminating reads. The used to run this analysis with paths to local `fastq` files is included at ['config/samples.csv`](config/samples.csv). Modifying that file with the changes below will run the analysis from the runs published on the SRA. 

To run this analysis with the data on the SRA, simply change the sample path in the [`config file`](/config/config.yml) from `config/samples.csv` to `config/samples_sra.csv`. 

| Run         | LibraryLayout | Virus | Host  | Source | Day |
|-------------|---------------|-------|-------|--------|-----|
| SRR13160722 | PAIRED        | SARS2 | human | public | 152 |
| SRR13160723 | PAIRED        | SARS2 | human | public | 146 |
| SRR13160724 | PAIRED        | SARS2 | human | public | 143 |
| SRR13160725 | PAIRED        | SARS2 | human | public | 130 |
| SRR13160726 | PAIRED        | SARS2 | human | public | 128 |
| SRR13160727 | PAIRED        | SARS2 | human | public | 81  |
| SRR13160728 | PAIRED        | SARS2 | human | public | 75  |
| SRR13160729 | PAIRED        | SARS2 | human | public | 25  |
| SRR13160730 | PAIRED        | SARS2 | human | public | 18  |

To run the analysis locally, you can use the following command. I would not recommend this, because the computational time will be extensive. 

```
snakemake --use-conda --conda-prefix ./env --cores 4
```

If you plan to run the analysis on Fred Hutch `rhino` or `gizmo` servers, you can submit the pipeline to `Slurm` with the following command. 

```
sbatch run_analysis.bash
```



