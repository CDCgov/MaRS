# MaRS (Malaria Resistance Surveillance)


The emergence of resistance to all currently available antimalarial drugs in multiple regions of the world represents a current global public health challenge. In order to monitor and address this situation, faster and more effective surveillance tools are required to track and monitor the emergence and evolution of drug resistance in malaria. The Malaria Resistance Surveillance (MaRS) project aims to address this challenge by collating and mapping genetic polymorphisms associated with drug resistance in malaria around the world. The project achieves this by employing a targeted amplicon deep sequencing (TADS) approach [Lab Protocol](https://github.com/CDCgov/MaRS/tree/master/lab_sop) to detect single nucleotide polymorphisms on all major malaria drug resistance genes associated genes in samples sourced from travelers returning to the US from overseas, as well as samples actively collected in collaboration with partners from other countries.

Data for this project can be found at the following link [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490). Collaborators are encouraged to submit their own data using this [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490)

The Malaria Resistance Surveillance or MaRS analysis pipline, is an attempt at standardizing the workflow for identifying both known and new polymorhisms in *P.falciparum* genes associated with drug resistance. 


## Version Histroy
  * Version 1.1.1 (1/9/2018)
    * Fastq samples detected from directory
    * Reference and dependencies for analysis on P.falciparum provided with the git bundle
    * BBDuk used to trim reads
    * BWA used for read alignment
    * Variant calling using Samtools and GATK
    * Summary tables for samples with separate tables for known and novel variants
    * Summary heatmaps and frequency graphs generated
    * The framework has been tested to work on Ubuntu 16.04 and RedHat Enterprise Linux Server Edition 6.8

# Setup

## Getting started

1. Download git repository:

Clone the master branch of this repository.
```{sh}
git clone https://github.com/CDCgov/MaRS.git
```

2. Download dependencies:

Download dependencies listed below. All the tools that are required to run MaRS included in the lib directory of the repository.
However, make sure that the Python and Java versions and packages are installed on the system.

3. Your first analysis:

Follow the directory structure listed below and use the run script included with the bundle to run your first analysis.
```{sh}
sh run.sh <path to experiment folder> <path to output folder>
```

For example if you have stored your fastq files in ```fq/``` folder and you want to store the results in the folder ```local/```. You can run the following command from the MaRS directory. 

```{sh}
sh run.sh fq/ local/
```

## Depdencies:

1. [Python3.4 ](https://www.python.org/download/releases/3.4.0/)
2. [Java (Version : 9.0.1)](http://download.oracle.com/otn-pub/java/jdk/9.0.1+11/jre-9.0.1_linux-x64_bin.tar.gz)
3. [Pandas (Version : 0.22.0)](http://pandas.pydata.org/pandas-docs/stable/)
4. [Numpy (Version : 1.13.3)](https://www.scipy.org/install.html)
5. [Seaborn (Version : 0.8.1)](https://seaborn.pydata.org/)
6. [Openpyxl (Version : 2.4.9)](https://pypi.python.org/pypi/openpyxl)
3. [BBMap (Version : v35.x)](https://sourceforge.net/projects/bbmap/)
4. [BWA (Version : 0.7.12)](http://bio-bwa.sourceforge.net/)
5. [Bowtie2 (Version : 2.3.3.1)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
6. [Snap (Version : 1.0beta23)](http://snap.cs.berkeley.edu/)
7. [Samtools (Version : 1.3.1)](http://www.htslib.org/)
8. [Bctools (Version : 1.3.1)](http://www.htslib.org/)
9. [GATK (Version : 3.6-0-g89b7209)](https://software.broadinstitute.org/gatk/download/)

## Directory structure:

1. fq

  * Contains all the input fastq files

2. lib

  * Contains the binaries for all the tools that MaRS can run

3. local

  * Local output folder

4. pyamd

  * Contains all the MaRS classes.

5. ref

  * Contains reference sequences. For the current version, the ref folder must also contain the alignment indicies.

## Creating alignment indicies:

1. BWA:
  ```{sh}
bwa index <reference-fasta.fa>
  ```


## Output directory structure:
* For each sample

  1. CleanedFastq: Folder containing adapter trimmed and quality filtered fastq files
  2. output.sam : Contains reads aligned to reference genome
  3. output_sorted.bam : Sorted BAM containing aligned reads
  4. output_sorted_RG.bam : Read gtroup added BAM file
  5. outout_fixmate.bam : Final BAM file, with mate information corrected
  6. variants.vcf : Variant calls from samtools
  7. variants_gatk.vcf : Variant calls from GATK HaplotypeCaller
  8. variants_samtools.bed : Annotated variant calls from samtools, in tab delimited format
  9. variants_gatk.bed : Annotated variant calls from GATK, in tab delimited format

* For the study

  1. Study_variants.xlsx : Summary table of all known variants that confer drug resistance, for all samples in the study
  2. Study_depth.xlsx : Summary of depth of coverage for codon correponding to variants that confer drug resistance
  3. Study_al_freq.xlsx : Summary of allele frequency of variants that confer drug resistance
  4. Study_novel_exonic_variants.xlsx : Summary of all novel variants found in exonic regions, for all samples in the study
  5. Study_novel_intronic_variants.xlsx : Summary of all novel variants found in intronic regions, for all samples in the study
  6. Study_novel_var_af.xlsx : Summary of allele frequency of novel variants for all samples
  7. Study_novel_var_depth.xlsx : Summary of depth of coverage for codon correponding to novel variants for all samples
  8. voi_af_heatmap.jpg : Heatmap showing the allele frequency across all samples for variants known to confer drug resistance
  9. voi_dp_heatmap.jpg : Heatmap showing the codon coverage across all samples for variants known to confer drug resistance
  10. voi_frequency.jpg : Count plot showing frequency variants conferring drug resistance across all the samples in the study
