


# MaRS (Malaria Resistance Surveillance)

The Malaria Resistance Surveillance or MaRS pipline, is an attempt at standardizing the workflow for identifying variants in parasite genomes conferring drug resistance. By streamlinging the data QC, alignment, and variant calling process, MaRS hopes to reduce time of analysis, as well allow for flexibilty of protocol used for the BGS analysis. The modular structure of the MaRS pipeline enables the user to build their own pipeline using the tools of their preference. The tools that MaRS is able to run are listed below, however, adding a custom tool the workflow is very simple. 

## Tools avaiable:

### Data quality analysis:
1. BBDuk

### Alignment methods:
1. BWA
2. Bowtie2
3. BBMap
4. Snap

### Variant calling:
1. Samtools
2. GATK

## Depdencies:

1. [Python3.4 ](https://www.python.org/download/releases/3.4.0/)
2. [BBMap](https://sourceforge.net/projects/bbmap/)
3. [BWA](http://bio-bwa.sourceforge.net/)
4. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
5. [Snap](http://snap.cs.berkeley.edu/)
6. [Samtools](http://www.htslib.org/)
7. [Bctools](http://www.htslib.org/)
8. [GATK](https://software.broadinstitute.org/gatk/download/)

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

1. Bowtie2:
  ```{sh}
bowtie2-build <refernce-fasta.fa> <reference-fasta.fa>
  ```

2. BWA:
  ```{sh}
bwa index <reference-fasta.fa>
  ```

3. Snap:
  ```{sh}
snap-aligner index <reference-fasta.fa> <output_directory>
  ```

## Output directory structure:
1. CleanedFastq: Folder containing adapter trimmed and quality filtered fastq files
2. output.sam : Contains reads aligned to reference genome
3. output_sorted.bam : Sorted BAM containing aligned reads
4. output_sorted_RG.bam : Read gtroup added BAM file
5. outout_fixmate.bam : Final BAM file, with mate information corrected
6. variants.vcf : Variant calls from samtools
7. variants_gatk.vcf : Variant calls from GATK HaplotypeCaller
8. variants_samtools.bed : Annotated variant calls from samtools, in tab delimited format
9. variants_gatk.bed : Annotated variant calls from GATK, in tab delimited format




