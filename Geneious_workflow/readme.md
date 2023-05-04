
<details>
  <summary><strong>Updates Over Time</strong></summary>

> Author: @ET 4/5/22 :goat:  
>> Edited: @MB 4/21/23 
----
>#### To Do ####

>#### Completed Activity ✓ ####

- [x] Updated readme to include Table of Contents
- [x] Updated readme to follow similar structure as all other readme files
- [x] Added MaRS Analysis workflow diagram to readme

------
</details>


# Overview of the Geneious_workflow

## Table of contents ## 


 * [Background](#Background)
 * [How do I get set up?](#Setup)
 * [Using the Geneious workflow](#GeneiousWorkflow)
 * [Contribution guidelines](#Guidelines)
 * [Who do I talk to?](#Contact)

 * * *


<a id="Background"></a>

## Background ##

The workflow described herein provides a step by step guide on how to analyze Illumina targeted deep amplicon (TADS) data for variants (SNPs) in drug resistance associated genes: _Pfk13_, _Pfdhfr_, _Pfdhps_, _Pfmdr1_, _Pfcrt_, and _Pfcytb_. Target audience are biologist or beginner bioinformaticans :beginner: :computer:


<a id="Setup"></a>

### How do I get set up? ###
* The user is encouraged to have a basic understanding of:
  * Scripting using [python](https://realpython.com/learning-paths/writing-pythonic-code/)
  * Working in a [terminal](https://mrkaluzny.com/blog/terminal-101-getting-started-with-terminal/) and [bash shell](https://linuxconfig.org/bash-scripting-tutorial-for-beginners)
  * Working with [jupyter-lab](https://pandas.pydata.org/getting_started.html) and [pandas](https://realpython.com/search?q=pandas)


* Required bioinformatics software:   
  * [Geneious Prime](https://www.geneious.com/prime/)


* Setting up your working environment:
  * [Virtual and python environment](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/py_ve_setup.md)
  * Additional required dependencies and modules will be listed in each jupyter-lab notebook
    * These you install in the terminal by `python3 -m pip install <module name>`


<a id="GeneiousWorkflow"></a>

## Using the Geneious workflow ##
The workflow is split into :four: parts:
* [01. Sample ID QC](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/)
* [02. Geneious Analysis](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis)
* [03. Summary Tables and Data Visualization](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_geneious_reports)
* [04. NCBI Submission](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/05_ncbi_submission)

The steps that make up the WorkFlow are outlined in **Figure 1**:

**Figure 1. Geneious Workflow Outline**

<img
  src="/images/MaRS_analysis_workflow.png"
  width="600">


<a id="Guidelines"></a>

### Contribution guidelines ###

* We highly encourage further improvment of this workflow  :thumbsup:
* Please fork :arrow_forward: make changes :arrow_forward: pull request
* Changes will be reviewed and tested before implementation
* We recommend setting up [git](https://www.atlassian.com/git) extension in jupyter-lab:
>   * Option 1: start `jupyter-lab`, click on jigsaw icon and enable extensions > search git > enable
>   * Option 2: `pip3 install --upgrade jupyterlab-git`, restart `jupyter-lab`



<a id="Contact"></a>

### Who do I talk to? ###
* We recommend posting questions or suggestions on
<a id="gitter"></a>
[![Gitter](https://badges.gitter.im/placeholder.svg)](https://gitter.im/placeholder?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
* We are a small group, however, we will do our best to respond in a timely fashion.  

:back: [To main page](https://github.com/CDCgov/MaRS)




