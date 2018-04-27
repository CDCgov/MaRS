# MaRS (Malaria Resistance Surveillance)


The emergence of resistance to all currently available antimalarial drugs in multiple regions of the world represents a current global public health challenge. In order to monitor and address this situation, faster and more effective surveillance tools are required to track and monitor the emergence and evolution of drug resistance in malaria. The Malaria Resistance Surveillance (MaRS) project aims to address this challenge by collating and mapping genetic polymorphisms associated with drug resistance in malaria around the world. The project achieves this by employing a targeted amplicon deep sequencing (TADS) approach [Lab Protocol](https://github.com/CDCgov/MaRS/tree/master/lab_sop) to detect single nucleotide polymorphisms on all major malaria drug resistance genes associated genes in samples sourced from travelers returning to the US from overseas, as well as samples actively collected in collaboration with partners from other countries.

Data for this project can be found at the following link [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490). Collaborators are encouraged to submit their own data using this [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490)

The Malaria Resistance Surveillance or MaRS analysis pipline, is an attempt at standardizing the workflow for identifying both known and new polymorhisms in *P.falciparum* genes associated with drug resistance.

*If you end up using MaRS in your workflow, please cite this [study](https://www.ncbi.nlm.nih.gov/pubmed/29439965):*

```
Next-Generation Sequencing and Bioinformatics Protocol for Malaria Drug Resistance Marker Surveillance.

Talundzic E1, Ravishankar S2, Kelley J3, Patel D4, Plucinski M4, Schmedes S4,5, Ljolje D3, Clemons B6, 
Madison-Antenucci S6, Arguin PM4, Lucchi NW4, Vannberg F2, Udhayakumar V4.

Antimicrob Agents Chemother. 2018 Mar 27;62(4). pii: e02474-17. doi: 10.1128/AAC.02474-17. Print 2018 Apr.
```

# Setup

## Getting started

1. Download git repository:

Clone the master branch of this repository.
```{sh}
git clone https://github.com/CDCgov/MaRS.git
```

2. Setup virtualenv for Python3:

MaRS requires python3 to be installed with pip available. Please make sure this is available on the system.
To avoid clashes with system version of required python modules, we recommend using a virtualenv
Run the following command to install virtualenv, if you already have virtualenv installed

```{sh}
python3 -m pip install virtualenv
virtualenv mars_env                   # Setup mars virtual environment
source mars_env/bin/activate          # Activate virtual environment
```
> If successfully activated, you should see now (mars_venv) in front of your terminal username.

3. Installing python dependencies and downloading third party libraries:

MaRS uses many python modules
Run the following command to install the dependencies
```{sh}
pip3 install pyvcf pysam matplotlib seaborn pandas numpy xlrd openpyxl
pip3 list --format=columns
```

4. Your first analysis:

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
3. [Pandas (Version : > 0.22.0)](http://pandas.pydata.org/pandas-docs/stable/)
4. [Numpy (Version : > 1.13.3)](https://www.scipy.org/install.html)
5. [Seaborn (Version : > 0.8.1)](https://seaborn.pydata.org/)
6. [Openpyxl (Version : > 2.4.9)](https://pypi.python.org/pypi/openpyxl)
3. [BBMap (Version : v35.92)](https://sourceforge.net/projects/bbmap/)
4. [BWA (Version : 0.7.12)](http://bio-bwa.sourceforge.net/)
5. [Bowtie2 (Version : 2.3.3.1)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
6. [Snap (Version : 1.0beta23)](http://snap.cs.berkeley.edu/)
7. [Samtools (Version : 1.3.1)](http://www.htslib.org/)
8. [Bcftools (Version : 1.3.1)](http://www.htslib.org/)
9. [GATK (Version : 3.6-0-g89b7209)](https://software.broadinstitute.org/gatk/download/)


## Version Histroy

  * Version 1.2.0 (2/28/2018)
    * Annotation script modified to account for MT genes.
  * Version 1.1.1 (1/9/2018)
    * Fastq samples detected from directory
    * Reference and dependencies for analysis on _P.falciparum_ provided with the Git bundle
    * BBDuk used to trim reads
    * BWA used for read alignment
    * Variant calling using Samtools and GATK
    * Summary tables for samples with separate tables for known and novel variants
    * Summary heatmaps and frequency graphs generated
    * The framework has been tested to work on Ubuntu 16.04 and RedHat Enterprise Linux Server Edition 6.8


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

## Output directory structure:
* For each sample

  1. CleanedFastq: Folder containing adapter trimmed and quality filtered Fastq files
  2. output.sam : Contains reads aligned to reference genome
  3. output_sorted.bam : Sorted BAM containing aligned reads
  4. output_sorted_RG.bam : Read group added BAM file
  5. outout_fixmate.bam : Final BAM file, with mate information corrected
  6. _Sample-name_\_variants_samtools_annotated.vcf : Annotated variant calls from Samtools
  7. _Sample-name_\_variants_gatk_annotated.vcf : Annotated variant calls from GATK HaplotypeCaller
  8. _Sample-name_\_variants_merge_annotated.vcf : Merged and annotated variant calls from GATK HaplotypeCaller and Samtools

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

## Public Domain
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This soruce code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.


## Privacy
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
Surveillance Platform [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
