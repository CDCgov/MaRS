# MaRS (Malaria Resistance Surveillance)

The emergence of resistance to all currently available antimalarial drugs in multiple regions of the world represents a current global public health challenge. In order to monitor and address this situation, faster and more effective surveillance tools are required to track and monitor the emergence and evolution of drug resistance in malaria. The Malaria Resistance Surveillance (MaRS) project aims to address this challenge by collating and mapping genetic polymorphisms associated with drug resistance in malaria around the world. The project achieves this by employing a targeted amplicon deep sequencing (TADS) approach [Lab Protocol](https://github.com/CDCgov/MaRS/tree/master/lab_sop) to detect single nucleotide polymorphisms on all major malaria drug resistance genes associated genes in samples sourced from travelers returning to the US from overseas, as well as samples actively collected in collaboration with partners from other countries.

Data for this project can be found at the following link [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490). Collaborators are encouraged to submit their own data using this [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490)

The NeST analysis framework was developed as the analysis framework for the MaRS project, in an attempt to standardizing the workflow for the identification both known and new polymorphisms in P.falciparum genes associated with drug resistance.

## Next-generation Sequence-analysis Toolkit (NeST) : A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data

Advancements in next-generation sequencing have led to the development of numerous bioinformatics tools and pipelines. Current tools for variant calling offer high-quality solutions; however, many tools are tailored for model organisms. Here, we present NeST, a consensus-based variant calling tool with a plug-and-play framework for use with any organism with minimal user input. NeST consists of four modules, integrating open-source bioinformatics tools and a custom VCF parser, to generate high-quality consensus variant calls. NeST was validated using targeted-amplicon deep sequencing data from 245 Plasmodium falciparum isolates to identify single-nucleotide polymorphisms conferring drug resistance. NeST offers a light-weight pipeline for variant calling with standardized outputs and minimal computational demands for easy deployment for use with various organisms and applications. The following document outlines details of installation, and usage of individual modules for analysis.

1. [Overview of NeST framework](#Overview)
2. [Availability of code and installation](#Installation)
3. [Input standardization](#inputs)
4. [NeST class structure](#classes)

<a id="Overview"></a>
### Overview of the NeST framework:

NeST is a python based modular framework for consensus based variant calling. The overall analysis framework is broken down into four major blocks.
1. PrepInputs
2. VarCallEngine
3. VCFToolkit
4. Summarize

![NeST framework overview](images/Kookaburra.png)

The figure outlines the four key blocks of NeST and the steps performed by each step. VarCallEngine and VCFToolkit are spawned in parallel for each sample that is being analyzed in the study. By default, 4 parallel threads are spawned, to account for minimum available computational resource. This can be altered as per availability of resources.

<a id="Installation"></a>
### Availability of code and installation:

1. Download git repository:

   Clone the master branch of this repository.
   ```
   git clone https://github.com/shashidhar22/NeST
   ```

2. Installing Python modules for NeST:
   NeST requires [Python3](https://www.python.org/downloads/) to be installed with [Pip](https://pip.pypa.io/en/stable/installing/) available. Please make sure this is available on the system. NeST uses many python modules, that can be installed using the following command

   ```
   python3 -m pip install numpy xlrd openpyxl pandas pysam --user
   ```

3. Installing R packages for Nest:

   NeST uses R for the visualization of the data, to enable this, the following modules need to be installed in R. 

   ```{R}
   install.packages(c('ggplot2', 'dplyr', 'tidyr', 'optparse', 'readr', 'stringr'))
   ```

4. Your first analysis:

   NeST comes packaged with an SRA accession list from the [MaRS]((https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490) experiment. The includes the SRA accession for 10 Illumina paired end samples. Running the command listed below, will download the 10 samples using [SRAToolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) and run NeST on it.

   ```
   sh run.sh fq/MaRS_test/SRR_Acc_List.txt local/MaRs_test
   ```

   To run NeST on locally stored fastq files. You can just provide the path to the input directory instead of the accession list.
   For example if you have stored your fastq files in ```fq/``` folder and you want to store the results in the folder ```local/```. You can run the following command from the NeST directory.

   ```
   sh run.sh fq/ local/
   ```

<a id="inputs"></a>
### Input standardization:

NeST is designed to reduce the amount of user intervention with regards to inputs that the user needs to provide. However to enable standardization of inputs across all organisms we require that a particular file format be followed for the three inputs listed below:

1. Fastq files:

   The PrepInputs module in NeST highly simplifies the management of fastq files. The module accepts two input formats.
   1. Input directory path:

      This just requires the user to provide the path to a folder containing fastq files. The files are recognized by the file extension, so the files must have either ```fq```, ```fq.gz```, ```fastq``` or ```fastq.gz``` file extensions. The name convention of paired file can be ```_1```, ```_r1```, or ```_R1```.

   2. SRA accession list:

      This list requires a ```.txt``` with a list of SRA experiments, with one SRA number per line. This can be export from the SRA run selector tool.
      An example SRA accession is provided under ```fq/MaRS_test/SRA_Acc_list.txt```.

2. BED format:

   The BED (Browser Extensible Data) is an easy and lightweight format to list annotations for a genome. NeST uses a full BED or BED 12 column format file as a guide to annotate variants with codon and amino acid changes. The example file listed below shows the details of how to structure the BED file. The separation of contig, gene and exon level information makes this format highly portable across genomes. The BED 12 column format for most organisms can be export from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). A detail explanation of the BED format can be found [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

   ```
   #contig start stop gene score strand  CDSstart  CDSstop rbg NoOfExons ExonLengths ExonStarts
   PfCRT 1	3399	PfCRT	.	+	95	3191	0	13	90,268,172,132,71,75,82,50,56,92,44,54,76,	96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115,
   MT	1	5967	COXIII	.	-	734	1573	0	1	839,	734,
   MT	1	5967	COL	.	+	1933	3471	0	1	1538,	1933,
   MT	1	5967	CYTOb	.	+	3492	4622	0	1	1130,	3492,
   PfDHFR	1	1827	PfDHFR	.	+	1	1827	0	1	1827,	1,
   PfDHPS	1	2417	PfDHPS	.	+	1	2417	0	3	135,1868,115,	1,313,2302,
   PfK13	1	2181	PfK13	.	+	1	2181	0	1	2181,	1,
   PfMDR1	1	4260	PfMDR1	.	+	1	4260	0	1	4260,	1,
   ```

3. Variant of Interest:

   The Summarize module within NeST, allows for easy summarization of variants called from all samples in a study. If a user specifies a list of variants of interest, a separate table will be created for these set of variants. The variants can be specified in ```.tsv```, ```.csv```, ```.xlsx``` format. And follows the format listed below

   | Chrom  | Gene   | RefAA | AAPos | AltAA |
   |:------:|:------:|:-----:|:-----:|:-----:|
   | PfCRT  | PfCRT  |   C   |   72  |   S   |
   | PfCRT  | PfCRT  |   V   |   73  |   V   |
   | PfMDR1 | PfMDR1 |   N   |   86  |   Y   |
   | PfMDR1 | PfMDR1 |   Y   |   184 |   F   |
   | MT     | CYTOb  |   I   |   258 |   M   |

## Citing MaRS 

*If you end up using MaRS in your workflow, please cite this [study](https://www.ncbi.nlm.nih.gov/pubmed/29439965):*

```
Next-Generation Sequencing and Bioinformatics Protocol for Malaria Drug Resistance Marker Surveillance.

Talundzic E, Ravishankar S, Kelley J, Patel D, Plucinski M, Schmedes S, Ljolje D, Clemons B,
Madison-Antenucci S, Arguin PM, Lucchi NW, Vannberg F, Udhayakumar V.

Antimicrob Agents Chemother. 2018 Mar 27;62(4). pii: e02474-17. doi: 10.1128/AAC.02474-17. Print 2018 Apr.
```
