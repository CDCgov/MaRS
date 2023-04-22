<details>
  <summary><strong>Updates Over Time</strong></summary>

> Author: @ET 4/5/22 :goat:  
>> Edited: @[MB]((mailto:mbajic@cdc.gov)) 4/21/23 
----
>#### To Do ####

>#### Completed Activity ✓ ####

- [x] Updated readme to include Table of Contents
- [x] Reformated readme to read easier with clear sections and steps
- [x] Specified two potential errors that can disrupt the workflow and how to address this
- [x] Added the current version of the MaRS Geneious workflow to the directory
- [x] Added the current version of the MaRS reference genes to the directory
- [x] Added the stable version of the Bowtie2 plugin to the directory for ease of access
- [x] Added a table listing the SNPs and haplotypes observed in commonly used _Plasmodium falciparum_ strains 
- [x] Added MaRS analysis workflow diagram to readme
- [x] Removed old Supplemental section and added new one that outlines what to do with the output from Geneious and how to evaluate the quality of the results

------
</details>



April 21st, 2023. Version 1.1.1; updated by [Marko Bajic](mailto:mbajic@cdc.gov)

> _Version 1.1: Updates made in the gene naming, delivery, and explanation of the Geneious workflow file "MaRS_Geneious_workflow_V1.geneiousWorkflow" and the reference files "MaRS_ReferenceGenes.geneious". Additionally, clarification about Bowtie2 version 2.3.0 was highlighted. The Supplemental section was updated to provide useful information for which step in the workflow to go to next, how to check the validity of positive and negative controls, and how to check the quality of the sequencing run using read numbers in the imported files, trimmed files, and alignment files.

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
* [Supplemental](#Supplemental)
    * [Next Steps](#NextgSteps)
    * [Evaluating Positive and Negative Controls](#Controls)
    * [Evaluating Quality of the Sequencing Run](#SeqQuality)
       * [Sequencing Quality Assessment Rationale](#SeqQRationale)
       * [Sequencing Quality Assessment Protocol](#SeqQProtocol)

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

The analysis steps relevant to Geneious are demonstrated in the yellow outline in **Figure 1**.

**Figure 1. MaRS Analysis Workflow**
<img
  src="/images/MaRS_analysis_workflow.png"
  width="600">


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

5. The MaRS_Geneious_workflow_V1 worklfow utilizes Bowtie2 for alignment. This aligner is not installed in Geneious by default. **Before the workflow can be run this plugin must be installed**. While this can be done through the Geneious Tools **some versions of Bowtie2 are not fully compatible with the MaRS workflow, or for aligning some files in Geneious in general**. For this reason, we provided the Bowtie2 version that works without issue with the workflow. This file is called [Bowtie_7_2_1.gplugin](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/Bowtie_7_2_1.gplugin) and is found in the [02_geneious_analysis](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/) folder; it was originally downloaded from [Geneious's plugins Download section](https://www.geneious.com/plugins/bowtie-plugin/#history). This plugin can be downloaded from GitHub to the user's computer. To install it in Geneious simply drag and drop the file into a Geneious window. A prompt will appear notifying that the plugin has been installed. Since the workflow is defined to use Bowtie2, as long as it is present as a plugin there is nothing that needs to be done to run the workflow. To specify, once the plugin has been installed by dragging and dropping it into Geneious the Workflow will be ready to run. The full list of all installed plugins can be seen in Geneious by navigating to Tools > Plugins. A notification from Geneious will pop up in the future with the notification that "A new plugin update is availble for download" but this should not be done. Instead click the "Don't check for new plugin updates" option and then choose "Don't Install".

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
    ✓Minimum Length: 100 bp

<a id="Aligning"></a>

### **B. Align/Assemble -> Map to Reference (Identify which reads belong to which amplicons and annotate them)** ###
- The workflow uses Map to Reference for algining reads to the reference
1. Dissolve contigs and re-assemble:<br />
Reference sequence: MaRS_ReferenceGenes (6 sequences) - MaRS<br />
    ✓ Assemble each sequence list separately
2. Method:<br />
    ✓ Mapper: Bowtie2<br />
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

# Supplemental #


<a id="NextgSteps"></a>

## Next Steps ##

After finishing the export of information from Geneious, the next step would be to navigate to [03_summary_reports](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_summary_reports) and follow the directions from there. This step would take the Geneious output and summarize and annotate it so that only the most essential information would be present for each SNP. This includes metadata information for the SNP based on which sample it originates from, but more importantly it identifies whether the SNP site is a mutation, to what degree (amount of reads that had the mutation), or whether the site is wild type and how much coverage was present at that site to identify this site as wild type. This is an important metric because not observing a mutation at a SNP site is not the same as observing that site as wild type; absence of reads does not mean that site is wild type, it just means that the site was not observed and cannot be identified as either mutant or wild type.


<a id="Controls"></a>

## Evaluating Positive and Negative Controls ##

There are a variety of [_Plasmodium falciparum_ control strains](https://www.beiresources.org/Catalog.aspx?f_instockflag=In+Stock&q=plasmodium+falciparum&page=1&f_displaysearchname=Nucleic%2bAcids) that can be ordered through [BEI](https://www.beiresources.org/) and used in the MaRS work. These strains will have specific mutations for each of the MaRS genes that are indicative of that strain. Typically, during the MaRS work, at least two different positive controls and one negative control are used per plate during PCR amplification, and therefore all downstream steps (amplicon combination, library preparation, library pooling, and sequencing). Whichever controls are chosen, they are replicated on different plates, but not on the same plate. The results of the replicates among the different plates for a given strain are expected to be the same. The aniticipated results for eight different _Plasmodium falciparum_ control strains are outlined in the file [SNPs_and_Haplotypes_ControlStrains.xlsx](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis/SNPs_and_Haplotypes_ControlStrains.xlsx). For negative controls, water is used instead of a DNA sample during PCR amplification of the amplicons. It is not uncommon to get some reads in these samples, but it is uncommon for the number to be above 100 reads. Furthermore, these reads are not expected to be of high quality and often cannot pass the BBDuk read trimming step because the output is empty. This is what is expected. Any reads that are of high quality are a source of contamination, but if the laboratory protocol is followed as described and special care is taken to prepare the PCR amplification mixtures in a separate hood from where DNA is added the potential for contamination is smaller than 1%. The expected amount of noise, low quality reads found in negative control samples, from all the negative control samples in a libray make up about 0.000006% of the library's reads on average. We will continue to monitor this number in the future to provide a more accurate estimation of this value.


<a id="SeqQuality"></a>

## Evaluating Quality of the Sequencing Run ##


<a id="SeqQRationale"></a>

### Sequencing Quality Assessment Rationale ###

Within Geneious, for any file that contains reads, the amount of reads present within the file can be identified in the column titled "# Sequences". This information is extremely important for evaluating the quality of the lab work as well as the sequencing run. The entire MaRS laboratory SOP can be evaluating by looking at the number of reads at these different stages of the Geneious analysis: 
 * Starting reads (RawReads: reads obtained for each sample during sequencing)
 * Reads after trimming (HQreads: high quality reads that remain after trimming by BBDuk)
 * Aligned reads (TotalAligned: amount of high quality reads that align to the MaRS genes). This can be further broken down to see how the distribution of each gene compares among the amplicons:
    * Reads that align to gene 1
    * Reads that align to gene 2
    * Reads that align to gene 3
    * Reads that align to gene 4
    * Reads that align to gene 5
    * Reads that align to gene 6

Overall, this information would look, for example, as shown in Table 1:

**Table 1. Read Metrics for Assessing Quality of Lab Work and Sequencing Run**

| AMDID                | RawReads | HQreads | PfdhpsReads | Pfk13Reads | Pfmdr1Reads | TotalAligned |
| ---                  | ---      | ---     | ---         | ---        | ---         | ---          |
| 21BFGO00A2015PfF3361 | 146932   | 86686   | 14619       | 38075      | 28452       | 20138        |
| 21BFGO00A2017PfF3361 | 130994   | 69916   | 13064       | 32726      | 20685       | 909          |
| 21BFGO00A2028PfF3361 | 131244   | 54454   | 21191       | 3727       | 26086       | 48188        |
| 21BFGO00A2034PfF3361 | 136766   | 70640   | 16240       | 35726      | 17120       | 56636        |
| 21BFGO00A2065PfF3361 | 109852   | 42170   | 10822       | 17472      | 12937       | 8821         |
| 21BFGO00A2066PfF3361 | 103812   | 54680   | 15471       | 23133      | 15088       | 39944        |
| 21BFGO00A2071PfF3361 | 145936   | 85380   | 21758       | 41076      | 20457       | 18257        |
| 21BFGO00A2087PfF3361 | 38166    | 12058   | 3304        | 8650       | 24          | 18138        |
| 21BFGO00A2089PfF3361 | 158818   | 44952   | 11768       | 22046      | 9927        | 77180        |
| 21BFGO00A2114PfF3361 | 144504   | 77286   | 20599       | 34546      | 20365       | 8430         |

 
Using the information from the RawReads column, the sum of these values for a sequencing run should give a number that is around the expected amount of reads afforded by the Illumina sequencing kits used for that library. For a V2 Regular kit, the expected amount of reads is 20 million, and for a V2 Nano kit the amount is 1 million reads. For the parameters set by the Geneious Workflow, there will be on average about 50% hiqh quality reads in the sequenced reads, and out of these around 97% will align to the MaRS genes that were amplified. As the sample data in Table 1 demonstrates, some MaRS genes will amplify with higher efficiency compared to other genes. In Table 1, these samples had about 45% of their high quality reads align to _Pfk13_ compared to 28% and 27% high quality reads that aligned to _Pfdhps_ and _Pfmdr1_, respectively. Further testing is needed to determine the overall amplification efficiencies for all six MaRS genes, but in general it is known that _Pfk13_ with better efficency. The amplification efficiency can be confirmed with gel electrophoresis and the GelAnalyzer software. If there is consistent amplification of one amplicon over the others, or the levels of amplification can be quanitified consistently among the different targets, then this can be used to determine at what level the amplicons should be combined before library purification: more amplified targets contribute a smaller volume during the combination of amplicons than targets that did not amplify as efficiently.


<a id="SeqQProtocol"></a>

### Sequencing Quality Assessment Protocol ###

To obtain the read numbers for each sample at the different stages of the Geneious workflow this process needs to be followed for all the files found in the three folders of interest, which are: 1) the starting folder where the files were imported, 2) the "Trimmed" folder where the BBDuk trimmed files are, and 3) the "Final_annotation" folder where the aligned files are located. The protocol for extracting the number of sequences (reads) in these files is:
1. At the top menu bar click File.
2. Navigate to Extract and click the "Documents..." option.
3. In the new prompt, select a location where to save the file to be exported.
4. Give the file to be exported a name, such as "{country name and year}RawReads.
5. For Files of Type choose "TSV tab-separated table (*.tsv).
6. Click Export.
7. In the new prompt, click "Proceed".
8. In the new prompt, select "Name" and "# Sequences".
9. Click "OK".

Currently, the exported information needs to be organized by hand using Excel, but we plan to write Python code to facilitate this proces in the future. 



