> [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)
>> Author: :baby_chick: JH  
>>> Edited: :goat: ET @ Date

> version 0.2

## Background ##

Molecular surveillance of anti-malarial resistance is an important public health activity. The workflow described herein provides a step by step guide on how to analyze Illumina targeted deep amplicon (TADS) data for variants (SNPs) in drug resistance associated genes: `k13, dhfr, dhps, mdr1, crt, and cytb`.

### How do I get set up? ###
* Required: [Geneious Prime](https://www.geneious.com/prime/)
* Basic scripting using [python](https://realpython.com/learning-paths/writing-pythonic-code/)
* Understanding of [jupyter-lab](https://pandas.pydata.org/getting_started.html) and [pandas](https://realpython.com/search?q=pandas)
* Required dependencies will be listed in each jupyter-lab notebook
* Test data is provided in [placeholder]()
* Target audience are biologist or beginner bioinformaticans :beginner: :computer:

### Contribution guidelines ###

* Further improvement of any of our workflows are _very_ encouraged :thumbsup:
* Please fork > make changes > pull request
* Changes will be reviewed and tested
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact



## ReadMe for Geneious workflow and scripts for reports for lab and EPI

####  1. Importing raw fastq files
        First import the raw fastq files to geneious. Then, select whole samples then click sequence on top menu. After clicking sequence, click set paired reads.
        For the option select paired reads and delete unpaired ones.
####  2. Geneious Workflow
    Geneious Workflow contains mainly four steps for processing raw sequences to variants and coverage. You select a samples then go to workflows and run the SNP_Coverage_Analysis.geneiousWorkflow. Then, it automatically runs all the below steps from a -> d.

###### a. Trimming (Discard unnecessary reads that would cause noise to the data)
    Geneious uses BBduk for trimming raw sequences. The settings are
    1. Trim Adapters: All Truseq, Nextera and PhilX adapters
    ✓ Trim Right End
    ✓ Kmer Length 27
    ✓ Maximum Substitutions 1
    ✓ Maximum Substitutions + INDELs: 0,
    2. Trim Low quality:
    ✓Both Ends
    ✓Minimum Quality 35
    3. Trim Adapters based on paired read overhangs: Minimum Overlap: 24
    4. Discard Short Reads: Minimum Length: 150 bp

###### b. Map to Reference (To observe differences in genetic structure from reference)
    The workflow uses Map to reference for algining reads to the reference
    1. Dissolve contigs and re-assemble: reference sequecne is updated_MaRS_refs_ET (6sequences),
    Assemble each sequence list separately
    2. Mapper:
    ✓ Bowtie2
    ✓ alginment Type: End to End
    ✓ Use Preset: High Sensitivity/Medium
    3. Do not trim (discrad trim annotations)
    4. Save contigs

###### c. Find Variations/SNPs (To find the unique genetic structure from the samples)
    1. Find Variants: Minimum Variant Frequency: 0.05
    2. Analtze effect of variants on translations : Default Genetic Code, Standard
    3. Advanced: Use separate annotations for each variant at a position
     ✓ Record names of all contributing sequences of each variant
     ✓ Only find variations in annotation types = TrackerSNP
     ✓ Also find varations within 0bases of those types
     ✓ CDS properties to Copy: gene, product, protein_id, locus, tag, note.

###### d. Find coverage (To see if there are WildTypes we are interested )
    1. Only Find in: Annotations in reference sequence of type: TrackerSNP, Create coverage annotations on reference sequence
    2. High Coverage: Find regions with coverage above, number of sequences: 0

###### e. Selecting workflow outputs (contigs )
    1. After the workflow is done contig documents will be generated in the same folder with AMD_ID and gene information
    2. To select the contigs you click on modified
    * If you can't see modified header scroll to the right
    3. Then select the files that are generated later with AMD_ID and gene information in the document title

###### f. Export findings to a csv raw output (To process more data )
    1. Select all the documents with coverage and variant information
    * Select in as mentioned in step e
    2. Select all the columns
    3. Select all the types
    4. Export table
    * Export as Annotations.csv


#### 2. Python scripts for Geneious raw input
###### a. Read Annotations.csv file and process
    1. Python script has with open function to open file and read each lines from the file
    2. Each jupyter note script will process Geneious raw output for different purposes.

###### b. File1: DF_Analysis1.ipynb/DF_Analysis1_EP.ipynb
    1. There are two versions one for lab and another one for EPI
    2. This version is written using pandas dataframe.
    3. The output gives information such as SITE,Sequence Name,G_Annotation,SNP,Mutation,WildType,VAF

####  3. test run  samples
        Naming schemia contains information  such as year, site, treatmentday, and pooled in order.
        First two number is year. Next two letter is country. Two letters after is site. The two numbers next means the treatment day.  
        The xp means it is pooled sample and not individual samples.
        List of example data for Angola with output from Geneious workflow.
        AN means angola, Be means Benguela, LS means Lunda Sul, Za means zaire
       Sample	Pooled	Year	SITE	TreatmentDay	GENE	G_annotation	COVERAGE	VAF	VF	SNP	Type	PooledSize	PooledGroup	PooledGroup
    1. 19ANBe00xp004PfF1152	pooled	19	Benguela	0	DHFR	S108N	1200	100.00%	1200	108N	mutation	10	4	4
    2. 19ANBe00xp005PfF1152	pooled	19	Benguela	0	DHFR	C59R	782	85.40%	668	59R	mutation	10	5	5
    3. 19ANBe00xp012PfF1152	pooled	19	Benguela	0	DHPS_437Corrected 	A437G	1072	100.00%	1072	437G	mutation	10	12	12
    4. 19ANBe00xp018PfF1152	pooled	19	Benguela	0	DHPS_437Corrected 	A437G	199	100.00%	199	437G	mutation	6	18	18
    5. 19ANLS00xp004PfF1152	pooled	19	Lunda Sul	0	DHPS_437Corrected 	S436A	1184	6.40%	76	436A	mutation	10	4	4
    6. 19ANLS00xp005PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	954	100.00%	954	108N	mutation	10	5	5
    7. 19ANLS00xp005PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	954	100.00%	954	108N	mutation	10	5	5
    8. 19ANLS00xp006PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	433	100.00%	433	108N	mutation	10	6	6
    9. 19ANZa00xp008PfF1152	pooled	19	Zaire	0	DHFR	C59R	485	96.50%	468	59R	mutation	10	8	8
    10. 19ANZa00xp009PfF1152	pooled	19	Zaire	0	PfCRT	I356T	122	15.60%	19	356T	mutation	10	9	9

    link to naming schema is: https://cdc.sharepoint.com/:p:/r/teams/CGH-DPDM-AMD/Shared%20Documents/SOPs%20-%20Lab/Sample%20Naming/Sample%20Naming%20Standardization_20211207.pptx?d=w0e0c338aa3a840e19a36e4684065345b&csf=1&web=1&e=pspDgd



    ### _Template Markdown_ ###

    ### What is this repository for? ###

    * Quick summary
    * Version
    * [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

    ### How do I get set up? ###

    * Summary of set up
    * Configuration
    * Dependencies
    * Database configuration
    * How to run tests
    * Deployment instructions

    ### Contribution guidelines ###

    * Writing tests
    * Code review
    * Other guidelines

    ### Who do I talk to? ###

    * Repo owner or admin
    * Other community or team contact
