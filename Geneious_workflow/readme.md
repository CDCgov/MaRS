## Background ##

The workflow described herein provides a step by step guide on how to analyze Illumina targeted deep amplicon (TADS) data for variants (SNPs) in drug resistance associated genes: `k13, dhfr, dhps, mdr1, crt, and cytb`. Target audience are biologist or beginner bioinformaticans :beginner: :computer:

### How do I get set up? ###
* You will need some basic understanding of:
  * Scripting using [python](https://realpython.com/learning-paths/writing-pythonic-code/)
  * Working in a [terminal](https://mrkaluzny.com/blog/terminal-101-getting-started-with-terminal/) and [bash shell](https://linuxconfig.org/bash-scripting-tutorial-for-beginners)
  * Working with [jupyter-lab](https://pandas.pydata.org/getting_started.html) and [pandas](https://realpython.com/search?q=pandas)


* Required bioinformatics software:   
  * [Geneious Prime](https://www.geneious.com/prime/)


* Setting up your working environment:
  * [Virtual and python environment](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/py_ve_setup.md)
  * Additional required dependencies and modules will be listed in each jupyter-lab notebook
    * These you install in the terminal by `python3 -m pip install <module name>`

### How do I get started ###
The workflow is split into :five: steps:
* [Sample ID QC](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/)
* [Geneious Analysis](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis)
* [Summary Tables](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_summary_tables)
* [Data Visualization](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/04_data_viz)
* [NCBI Submission](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/05_ncbi_submission)

### Contribution guidelines ###

* We highly encourage further improvment of this workflow  :thumbsup:
* Please fork :arrow_forward: make changes :arrow_forward: pull request
* Changes will be reviewed and tested before implementation
* We recommend setting up [git](https://www.atlassian.com/git) extension in jupyter-lab:
>   * Option 1: start `jupyter-lab`, click on jigsaw icon and enable extensions > search git > enable
>   * Option 2: `pip3 install --upgrade jupyterlab-git`, restart `jupyter-lab`


### Who do I talk to? ###
* We recommend posting questions or suggestions on
<a id="gitter"></a>
[![Gitter](https://badges.gitter.im/placeholder.svg)](https://gitter.im/placeholder?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
* We are a small group, however, we will do our best to respond in a timely fashion.  

:back: [To main page](https://github.com/CDCgov/MaRS)



----
<details>
  <summary><strong>TODO</strong></summary>

> Author: @ET 4/5/22 :goat:  
>> Edited: @ET 11/30/22
----
>#### TODO ####
>#### Activity Name ####

>#### Completed Activity ✓ ####

- [x] Add readme to each of the steps @ET :goat:
- [x] Add gitter link @ET :goat:  
- [x] Add links to each of the steps @ET :goat:
- [x] Update README: add python, pyenv, ven set up  @ET :goat:

------
</details>

