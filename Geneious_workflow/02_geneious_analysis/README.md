April 19th, 2023. Version 1.1.0; updated by [Marko Bajic](mailto:mbajic@cdc.gov)

> _Version 1.1: Updates made in the gene naming, delivery, and explanation of the Geneious workflow file "MaRS_Geneious_workflow_V1.geneiousWorkflow" and the reference files "MaRS_ReferenceGenes.geneious"

  * * *
# Analyzing NGS Sequences with Geneious Using the MaRS Workflow and Annotated Reference Files

## Table of contents ## 


* [Introduction](#intro)
    * [Background](#Background)
    * [Set-up](#Set-up)
    * [Contribution guidelines](#Contribution_guidelines)
    * [Who do I talk to?](#Who_do_I_talk_to)
* [Protocol](#Protocol)
    * [0. Utilizing test data to test the workflow](#UtilizingTestData)
    * [I. Importing raw fastq files](#ImportingFastq)
    * [II. Geneious Workflow](#GeneiousWorkflow)
    * [III. Procedural Steps for the MaRS Geneious Workflow](#Procedure)
        * [A. Trim using BBDuk](#Trimming)
        * [B. Align/Assemble -> Map to Reference](#Aligning)
        * [C. Find Variations/SNPs](#FindingSNPs)
        * [D. Find High/Low Coverage](#FindingCoverage)
        * [E. Export findings to a csv raw output](#ExportingGeneiousResults)
* [Supplemental: Python scripts for Geneious raw input](#Supplemental)


 * * *

<a id="intro"></a>

# Introduction #

<a id="Background"></a>
## Background ##

This is the standard Operating Procedure (SOP) for using the Geneious application to analyze Next-Generation Sequencing (NGS) results generated in lab using the [MaRS Laboratory SOP](https://cdcgov.github.io/MaRS/) for six full length _P. falciparum_ genes associated with antimalarial resistance:
* [_Plasmodium falciparum_ kelch 13 (_Pfk13_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_1343700)
* [_Plasmodium falciparum_ chloroquine resistant transporter (_Pfcrt_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0709000)
* [_Plasmodium falciparum_ multidrug resistant protein 1 (_Pfmdr1_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0523000)
* [_Plasmodium falciparum_ bifunctional dihydrofolate reductase-thymidylate synthase (_Pfdhfr_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0417200)
* [_Plasmodium falciparum_ dihydropteroate synthase (_Pfdhps_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0810800)
* [_Plasmodium falciparum_ cytochrome b (_Pfcytb_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_MIT02300)

The instructions outlined below are for setting up the two required files, the MaRS Geneious workflow and the reference file that contains the six genes analyzed commonly by the MaRS protocol, and for using the two files to analyze, from beginnig to the end, NGS sequences using Geneious. The results can then be exported for downstream analysis and reporting. Of particular note, the reference file contains annotations for each SNP of interest that is commonly tracked withint the MaRS lab. The specific track that annotates these SNPs of interest is named, within Geneious, the TrackerSNP. 

The information presented here is organized into steps and associated files contained within [02_geneious_analysis](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis). These steps follow after [01_sampleID_QC](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/01_sample_ID_QC) and precede [03_summary_tatbles](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_summary_tables). All of these procedural steps are organized within [Geneious_workflow](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow).



<a id="Set-up"></a>
## Set-up ##

* Required: [Geneious Prime](https://www.geneious.com/prime/)
* Basic scripting using [python](https://realpython.com/learning-paths/writing-pythonic-code/)
* Understanding of [jupyter-lab](https://pandas.pydata.org/getting_started.html) and [pandas](https://realpython.com/search?q=pandas)
* Required dependencies will be listed in each jupyter-lab notebook
* Test data will need to be downloaded using the [data_download jupyter notebook](https://github.com/CDCgov/MaRS/blob/goat_dev/Geneious_workflow/Test_data/data_download.ipynb) in the [Test data directory](https://github.com/CDCgov/MaRS/tree/goat_dev/Geneious_workflow/Test_data)  
* Target audience are biologist or beginner bioinformaticans

<a id="Contribution_guidelines"></a>
## Contribution guidelines ##

* Further improvement of any of our workflows are _very_ encouraged :thumbsup:
* Please fork > make changes > pull request
* Changes will be reviewed and tested
* Other guidelines

<a id="Who_do_I_talk_to"></a>
## Who do I talk to? ##

* Repo owner or admin
* Other community or team contact

* * *

<a id="Protocol"></a>

# Protocol #


<a id="UtilizingTestData"></a>

## 0. Utilizing test data to test the workflow

If you plan to use test data, please download the [Test data directory](https://github.com/CDCgov/MaRS/tree/goat_dev/Geneious_workflow/Test_data), start up [jupyter-lab](https://pandas.pydata.org/getting_started.html) and [pandas](https://realpython.com/search?q=pandas) and follow the directions in the [data_download jupyter notebook](https://github.com/CDCgov/MaRS/blob/goat_dev/Geneious_workflow/Test_data/data_download.ipynb).   


<a id="ImportingFastq"></a>

##  I. Importing raw fastq files ##

1. First import the raw fastq files to Geneious either by dragging and dropping the files from your computer into a Geneious folder or by selecting File > Import > Files from the menu.

2. Select all samples.

3. Click sequence on top menu; after clicking sequence, click set paired reads.

4. Under option select paired reads and delete unpaired ones.


<a id="GeneiousWorkflow"></a>

##  II. Geneious Workflow ##
- The Geneious Workflow ["MaRS_Geneious_workflow_V1"](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/MaRS_Geneious_workflow_V1.geneiousWorkflo) contains mainly four steps for calling SNPs. 
    1. Trimming of short and/or low quality reads.
    2. Aligning trimmed reads to reference genes.  
    3. Finding Variations/SNPs in aligned files at TrackerSNP sites.  
    4. Finding High/Low Coverage at aligned TrackerSNP sites.  

1. To add the [MaRS_Geneious_workflow_V1](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/MaRS_Geneious_workflow_V1.geneiousWorkflo) workflow to Geneious:
    - Download it from Github using the link and right clicking the "Raw" button and choosing "Save Link As..." option and a location on your computer and clicking "Save".
    - Find the file on your computer.
    - Drag and drop it into the Geneious application.

2. The MaRS_Geneious_workflow_V1 worklfow can be viewed and edited within Geneious by navigating to it via "Workflows" > "Manage Workflows...". 

3. In the workflow, the step "Align/Assemble -> Map to Reference(s)" specifies which reference genes to align trimmed reads to. By default, this is set to "MaRS_ReferenceGenes". It is is important to either have this reference file within the same folder where the Workflow is being run or to specify its location within Geneious by double clicking this step in the Workflow and changing the "Reference Sequence:" file by clicking "Choose...". 

4. The reference file [MaRS_ReferenceGenes.geneious](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/MaRS_ReferenceGenes.geneious) which contains the sequences, annotations, and TrackerSNPs for the six MaRS genes: _Pfk13_, _Pfdhfr_, _Pfdhps_, _Pfmdr1_, _Pfcrt_, and _Pfcytb_ is found on Github and can be downloaded using the same steps as outlined above for downloading the Workflow.
    - Download it from Github using the link and right clicking the "Raw" button and choosing "Save Link As..." option and a location on your computer and clicking "Save".
    - Find the file on your computer.
    - Drag and drop it into the Geneious application into a specific folder.
        - Alternatively, you can import it into Geneious trhough File > Import > Documents

5. The MaRS_Geneious_workflow_V1 worklfow utilizes Bowtie2 for alignment. This aligner is not installed in Geneious by default. Before the workflow can be run this plugin must be installed. This can be done by navigating through this menu in Geneious: Tools > Plugins. Click the "Install" button on the right side of the Bowtie option under "Available Plugins". Once done, the option "Bowtie short read mapper: Run bowtie plugin" should appear under the "Installed Plugins" section. Once this is done, the workflow can be run without issue.

6. To run the MaRS_Geneious_workflow_V1 worklfow, select the fastq files you would like to analyze and choose "Workflows" > "MaRS_Geneious_workflow_V1"
    - A prompt will ask for the maximum memory to use
        - Choose an appropriate value based on the amount of RAM available on the computer and how many other processes need to be run besides Geneious.
    - After pressing "OK", the workflow will run all the steps (from A to D) as outlined in the "III. Procedural Steps of the MaRS Geneious Workflow" section.

<a id="Procedure"></a>

##  III. Procedural Steps for the MaRS Geneious Workflow ##

<a id="Trimming"></a>

### **A. Trim using BBDuk (Discard low quality reads that would create noise and uncertainty in the results)** ###
- Geneious uses BBduk for trimming raw sequences. The settings are
- **NOTE:** This step is the only step that is known to be able to stop the workflow and produce no outputs. This occurs if at least one of the selected files to be analyzed by the worfklow contains 100% low-quality reads (no reads with minimum quality of 35 or higher and mininmum lenght of 150 bp or more). If this occurs and the "BBDuk produced no results" error is encountered there is an approach for how to identify which files are creating the error. The most likely candidate files that caused the error are those with low read amounts, such as a file with 200 reads or fewer. This information can be seen for each file in the "# Sequences" column. Excluding these samples will allow the workflow to be run without issue. Alternatively, each file can be run individually through the workflow and for any file that fails to run the workflow must be excluded as that file has insufficient reads of high quality.

1. Trim Adapters: All Truseq, Nextera and PhilX adapters<br />
    ✓ Trim Right End<br />
    ✓ Kmer Length 27<br />
    ✓ Maximum Substitutions 1<br />
    ✓ Maximum Substitutions + INDELs: 0<br />
2. Trim Low quality:<br />
    ✓Both Ends<br />
    ✓Minimum Quality: 35<br />
3. Trim Adapters based on paired read overhangs: <br />
    ✓Minimum Overlap: 24
4. Discard Short Reads:<br />
    ✓Minimum Length: 150 bp

<a id="Aligning"></a>

### **B. Align/Assemble -> Map to Reference (Identify which reads belong to which amplicons and annotate them)** ###
- The workflow uses Map to Reference for algining reads to the reference
1. Dissolve contigs and re-assemble:<br />
Reference sequence: MaRS_ReferenceGenes (6 sequences) - MaRS<br />
    ✓ Assemble each sequence list separately
2. Method:<br />
    ✓ Mapper: Geneious<br />
    ✓ Alignment Type: End to End<br />
    ✓ Use Preset: High Sensitivity/Medium
3. Do not trim (discard trim annotations)
4. Save contigs

<a id="FindingSNPs"></a>

### **C. Find Variations/SNPs (Compared to the reference, identify sequence variations in the assembled contig)** ###
 1. Find Variants: <br />
 ✓ Minimum Variant Frequency: 0.05
 2. Analtye effect of variants on translations:<br />
 ✓ Default Genetic Code, Standard
 3. Advanced: Use separate annotations for each variant at a position<br />
 ✓ Record names of all contributing sequences of each variant<br />
 ✓ **Only find variations in annotation types: TrackerSNP**<br />
 ✓ Also find varations within 0 bases of those types<br />
 ✓ CDS properties to Copy: gene, product, protein_id, locus, tag, note.

<a id="FindingCoverage"></a>

### **D. Find High/Low Coverage (Determine if the absence of a SNP is due to no coverage or because the site is Wild Type)** ###
1. Only Find in:<br />
 ✓ Annotations in reference sequence of type: TrackerSNP<br />
 ✓ Create coverage annotations on reference sequence
2. High Coverage:<br />
 ✓ Find regions with coverage above<br />
 ✓ number of sequences: 0


<a id="ExportingGeneiousResults"></a>

### **E. Export findings to a csv raw output (Use output to summarize essential results and generate a report)** ###
1. The alignment files with called SNPs and coverage information will be saved in a sub-folder called "Final_annotation".
2. Inside this sub-folder, select all the documents and in the bottom of the Geneious window select "View Documents"
    * **NOTE:** The "View Documents" option will only appear if a moderate amount of samples are chosen. If only a dozen samples are selected this prompt will not appear. In this case, just proceed to the next step.
3. In the window at the bottom of Geneious, follow these steps:
     1. Under Types, select Show one > Polymorphism.
     2. Then again, under Types choose "Coverage - High".
        * This should result in both "Polymorphism" and "Coverage - High" being selected for Types. You can visually confirm this by clicking Types and verifying that only these two options have a checkmark next to them.
     3. Under Columns select "Show All".
        * The "Show All" is the easiest and quickest way to select all of the required information, but it is more information than is needed and can theoretically cause some issues with some Jupyter setups on Windows.
        * If there is time to be diligent, instead of selecting "Show All", these twelve columns are the essential information from Columns that needs to be exported and can be selected one by one: Document Name, Sequence Name, Type, Amino Acid Change, Average Coverage, CDS Codon Number, CDS Position Within Codon, Coverage, gene, TrackerSNP, Variant Frequency, Variant Raw Frequency.
     4. Click "Export table".
     5. Select a location where to save the file and click "Save".
        * **Under "Files of Type:" choose "Comma Separated Values (*.csv)".**

 * * *

<a id="Supplemental"></a>
# Supplemental: Python scripts for Geneious raw input

## A. Read Annotations.csv file and process
1. Python script has with open function to open file and read each lines from the file
2. Each jupyter note script will process Geneious raw output for different purposes.

## B. File1: DF_Analysis1.ipynb/DF_Analysis1_EP.ipynb
1. There are two versions one for lab and another one for EPI
2. This version is written using pandas dataframe.
3. The output gives information such as SITE,Sequence Name,G_Annotation,SNP,Mutation,WildType,VAF

##  C. Test run  samples
1. Naming schemia contains information  such as year, site, treatmentday, and pooled in order.
2. First two number is year. Next two letter is country. Two letters after is site. The two numbers next means the treatment day.  
3. The xp means it is pooled sample and not individual samples.
4. List of example data for Angola with output from Geneious workflow.
5. AN means angola, Be means Benguela, LS means Lunda Sul, Za means zaire

       
        Sample	Pooled	Year	SITE	TreatmentDay	GENE	G_annotation	COVERAGE	VAF	VF	SNP	Type	PooledSize	PooledGroup	PooledGroup
        1. 19ANBe00xp004PfF1152	pooled	19	Benguela	0	DHFR	S108N	1200	    100.00%	1200	108N	mutation	10	4	4
        2. 19ANBe00xp005PfF1152	pooled	19	Benguela	0	DHFR	C59R	782	85.    40%	668	59R	mutation	10	5	5
        3. 19ANBe00xp012PfF1152	pooled	19	Benguela	0	DHPS_437Corrected 	    A437G	1072	100.00%	1072	437G	mutation	10	12	12
        4. 19ANBe00xp018PfF1152	pooled	19	Benguela	0	DHPS_437Corrected 	    A437G	199	100.00%	199	437G	mutation	6	18	18
        5. 19ANLS00xp004PfF1152	pooled	19	Lunda Sul	0	DHPS_437Corrected 	    S436A	1184	6.40%	76	436A	mutation	10	4	4
        6. 19ANLS00xp005PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	954	100.    00%	954	108N	mutation	10	5	5
        7. 19ANLS00xp005PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	954	100.    00%	954	108N	mutation	10	5	5
        8. 19ANLS00xp006PfF1152	pooled	19	Lunda Sul	0	DHFR	S108N	433	100.    00%	433	108N	mutation	10	6	6
        9. 19ANZa00xp008PfF1152	pooled	19	Zaire	0	DHFR	C59R	485	96.50%	    468	59R	mutation	10	8	8
        10. 19ANZa00xp009PfF1152	pooled	19	Zaire	0	PfCRT	I356T	122	15.    60%	19	356T	mutation	10	9	9

6. The link to the naming schema is: https://cdc.sharepoint.com/:p:/r/teams/CGH-DPDM-AMD/Shared%20Documents/SOPs%20-%20Lab/Sample%20Naming/Sample%20Naming%20Standardization_20211207.pptx?d=w0e0c338aa3a840e19a36e4684065345b&csf=1&web=1&e=pspDgd


----
<details>
  <summary><strong>TODO</strong></summary>


> Author: :baby_chick: JH @ 04/7/22 
>> Edited & Reviewed: :tiger: DP & :goat: @ET 11/30/22 
----
>#### TODO ####
>#### Activity Name ####

  - [ ] template action @XX  

>#### Completed Activity ✓ ####
  - [x] Update readme to point to new example fastq (external) @ET
   - [x] Add directions on how to download from SRA @DP

------
</details>
