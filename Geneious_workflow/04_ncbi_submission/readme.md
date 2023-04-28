
<details>
  <summary><strong>Updates Over Time</strong></summary>

> Author: @ET 4/5/22 :goat:  
>> Edited: @[MB]((mailto:mbajic@cdc.gov)) 4/21/23 
----
>#### To Do ####

>#### Completed Activity ✓ ####

- [x] Created readme file
- [x] Added two example files to the directory for demonstration of file organization relevant to NCBI submission

------
</details>


# Submitting Sequence Files to an Already Existing Bioproject #

## Table of contents ##  


* [Outline](#Outline)
  * [Rationale](#Rationale)
  * [Needed items](#NeededItems)
  * [Protocol Outline](#ProtocolOutline)
* [Protocol](#Protocol)
  * [I. Preparation](#Preparation)
  * [II. Submitting a new BioSample](#SubmitAnewBioSample)
  * [III. Submitting a new Sequence Read Archive (SRA)](#SubmitAnewSRA)
  * [IV. Confirming Completion](#CompletionSteps)
* [Supplemental](#Supplemental)
  * [Example Documents](#ExampleDocs)

 * * *


<a id="Outline"></a>

# Outline #

<a id="Rationale"></a>

## Rationale ##
The Malaria Resistance Surveillance (MaRS) project aims to address this challenge by collating and mapping genetic polymorphisms associated with drug resistance in malaria around the world. As samples are collected and sequenced, they need to be shared for public access, but more importantly for efficient sharing of sequencing data with our collaborators and team members. To that end, we created on NCBI the bioproject [PRJNA428490](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490) as the central NCBI repository for all of the sequencing data that is collected and shared by the MaRS team. Over time, new domestic and international samples will be sequenced and shared through this bioproject, instead of being shared as separated bioprojects. Therefore, it is important to know the steps involved in submitting newly sequenced files, and their accompanying metadata, to an already established bioproject, and that is what the guidelines below will describe.

<a id="NeededItems"></a>

## Needed Items ##
Ensure that these items are ready before proceeding to the using this protocol.
* NCBI account 
* Sequence read files, most commonly exist as “.fastq.gz” files
* Metadata sheet describing each sample (example sheet). This will be submitted as a tsv file but can be made in Excel and then exported as tsv (tab separated value) file.

<a id="ProtocolOutline"></a>

## Protocol Outline ##
0.  Preparation
1.	Submit a new BioSample
2.	Organize the read files
3.	Create a metadata sheet and fill out most of the information
4.	Submit a new Sequence Read Archive submission


<a id="Protocol"></a>

# Protocol #


<a id="Preparation"></a>

## I. Preparation ##
1. Organize your fastq files on your computer for upload to NCBI.
2. Fill out the metadata sheets for all of your fastq files. This may take you some time so finish this section before proceeding to the next section. **Please note, there will be TWO metadata files you will submit overall. One metadata file will be for the BioSample and the other will be for the SRAs.**
    1. For MaRS work, the template for the Biosample metadata sheet will be the one for microbes and it has several columns that are mandatory, some that only one value is mandatory, and many optional options. 
    
    2.	When submitting SRAs you will also have to provide a metadata sheet for each fastq file. However, you cannot fill out ALL of the information for the SRA metadata until AFTER you have submitted the Biosample information as one of the pieces of information you will have to provide will ask for the BioSample Accession that you will not have until after the BioSample information is submitted. Also, please note, paired-end reads will have information for 2 columns “filename” and “filename2” for each sample’s “sample_name” column. More on this in section III-5 below. A sample SRA Metadata sheet can be found here.


<a id="SubmitAnewBioSample"></a>

## II. Submitting a new BioSample
1.	Go to your [Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/),
2.	Choose “BioSample” from “Start a new submission”.
3.	Under **<ins>SUBMITTER</ins>** fill out the relevant information. Choose “No group.” For Submitting organization, Department, and Address use the appropriate inforamtion for your affiliation.
4.	Under **<ins>GENERAL INFO</ins>** for most projects you will choose “Release immediately following processing.” When submitting multiple samples you will choose Batch/Multiple BioSamples. Please note, this is where you can find the link to the templates page. 
5.	Under **<ins>SAMPLE TYPE</ins>** input “Plasmodium falciparum” or appropriate species, and then click “show all 20 packages” and choose Microbe if it is a Plasmodium sample.
6.	Under **<ins>ATTRIBUTES</ins>** choose “Upload a file using Excel or text format…” and then upload the corresponding Metadata sheet you have filled out for the Biosample.
7.	Under **<ins>REVIEW & SUBMIT</ins>** check to make sure everything looks correct, and press Submit.
8.	Await an email confirming processing has finished and this information is ready.


<a id="SubmitAnewSRA"></a>

## III.	Submitting a new Sequence Read Archive (SRA) ##
1.	Go to your [Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/).
2.	Choose “Sequence Read Archive” from “Start a new submission”.
3.	Under **<ins>SUBMITTER</ins>** fill out the relevant information, as in previous section (II-3).
4.	Under **<ins>GENERAL INFO</ins>** choose Yes for the “Did you already register a BioProject for this research…” question. This is IMPORTANT! For the “Existing BioProject” input “PRJNA428490” (no quotes). Under “Did you already register a BioSample for this sample…” choose Yes. Under “When should this submission be released to the public” choose “Release immediately following processing.”
5.	Under **<ins>SRA METADATA</ins>** you will choose “Upload a file using Excel or text format” from “How do you want to provide your metadata” and then upload the corresponding Metadata sheet for SRA. This was briefly addressed in B-2, but at this point you will have all of the information needed to fill out the metadata sheet, notably the BioSample Accession information. 

    - A.	To actually obtain this information, navigate to your [“My submission” page](https://submit.ncbi.nlm.nih.gov/subs/) and click the appropriate BioSample submission (leftmost option, under Submission). On the new page you will be able to click the “Download attributes file with BioSample accessions” and get the accession ID for each sample form there. Alternatively, you can go to Manage Data > BioSample > scroll to bottom and click Download [] records, but please note that this will include all your biosamples over time.
    - B.	You HAVE to name your sheet (if Excel) as “SRA_data” when submitting SRA metadata as an Excel file. If you do not, an error will prevent you from uploading. To that end, it is important that the headers are EXACTLY as NCBI wants them. This includes case and the use of “_” instead of a space in the headers. 
6.	Under **<ins>FILES</ins>** you will upload your actual fastq files, but this can be finnicky and really depends on how many (overall size) of your fastq files. It is the easiest option. Most often, you will submit files over the “FTP or Aspera Command Line file preload” option. Select this and click the “Request preload folder” option, it may ask you whether you’re sure you want to leave the page, just press Leave because it will bring you right back to the same page but now it will have instructions for “Aspera command line upload instructions” or “FTP upload instructions.” Command line (terminal) can be confusing, and so that windows users don’t have to jump too many hurdles the simplest path (Graphical User Interface of an FTP client) will be described:
    - A.	Go to https://filezilla-project.org/ and download the appropriate OS’s version of FileZilla (for me it was Windows). This is more user friendly than CyberDuck because of NCBI’s rather CONFUSING approach for connecting, as noted in “E” below.
    - B.	Install FileZilla
    - C.	Launch/open FileZilla
    - D.	At the top under Host, Username, and Password, provide the information as specified by NCBI on the FILES page we are still working on. Note, Host=Address
    - E.	Press Quickconnect and if you get "550 /: Permission denied" or "Failed to read the directory listing" just know that you did everything correctly. Weird indeed.
    - F.	On the right side of the window, under Remote site type what it says under step 3 from NCBI, for me it was “uploads/marko.bajic25_gmail.com_qzu1xigx” and press Enter. Now, on the left-hand side (Local site) navigate to where the fastq files are. Once that folder is opened the files will be visible in the lower-most left panel.
    - G.	THIS IS IMPORTANT! Before transferring files, right click either on your root folder (for me it was marko.bajic25_gmail.com_qzu1xigx) in the right middle pane, or the right bottom pane. Create a new folder and name it something appropriate. Then you can drag and drop all the fastq files into that folder. Do NOT drag and drop them into the root folder (I mean, technically, you can and then you create a subfolder and drag them in there from the root folder). Just make sure that when finished with these steps the fastq files are in a subfolder on NCBI’s end.
    - H.	To KNOW that this has worked, you need to give it 10 minutes. Then on the NCBI page (we are still on FILES part of SRA submission) the “Select preload folder” option will (after 10 minutes since FTP upload) have (you can press “refresh folders” if it does not) the folder you created in FileZilla listed. Now just select it and choose “Use selected folder.” PRO TIP: check the size of the folder you created and the size of that same folder on your local computer. If there is a big discrepancy then you may want to redo the submission of files, but if they are generally close (2.9 and 2.93 for me) then this acts as a sanity check. Same way, you can check that the number of files matches on both ends. 
    - I.	Scroll to bottom of page (NCBI, FILES) and click Continue
7.	It will take some time for the **<ins>FILES</ins>** step to finish PROCESSING. When it is finished, on the **<ins>REVIEW & SUBMIT</ins>** page you can visually inspect everything. Most notably, does each biosample accession have two files in paired end sequencing. When done click Submit. 
8.	Await email confirming processing is finished.


<a id="CompletionSteps"></a>

## IV.	Confirming Completion ##
1. Go to the MaRS Bioproject [PRJNA428490](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490).
2. Confirm that the new samples are there, and they are accurately annotated.


 * * *

<a id="Supplemental"></a>

# Supplemental #


<a id="ExampleDocs"></a>

## Example Documents ##
For convenience, two example files, [Example_Biosample_metadata.xlsx](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/Example_Biosample_metadata.xlsx) and [Example_SRAmetadata.xlsx](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/Example_SRAmetadata.xlsx), have been included in the [04_ncbi_submission](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/04_ncbi_submission) folder. These files contain 10 rows of information that gets submitted to NCBI. One file demonstrates how the Biosample metadata sheet can be filled out and the other demonstrates how the SRA metadata sheet can be filled out.


