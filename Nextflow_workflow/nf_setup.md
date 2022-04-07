> Author: @ ET updated 4/4/22 :goat:  
> version 0.2

## Background ##

The [nf-NeST](https://github.com/CDCgov/Nf-NeST) workflow uses [nextflow](https://www.nextflow.io/docs/latest/index.html) to analyze Illumina targeted deep amplicon (TADS) data for variants (SNPs) in drug resistance associated genes: `k13, dhfr, dhps, mdr1, crt, and cytb`. Target audience are advanced Bioinformaticans. :octocat::computer:

## What is this repo for? ##

* Help with setting up [nextflow](https://www.nextflow.io/docs/latest/index.html) to run the automated [nf-NeST](https://github.com/CDCgov/Nf-NeST) pipeline.

#### Nextflow ####

[Nextflow](https://www.nextflow.io/docs/latest/index.html) is a bioinformatics workflow management software that enables scalable and reproducible workflows. It works with Docker, Singularity and Conda. It uses a common scripting language and processes, allowing for easy integration and running on cloud or on-premise HPC environments. [Nf-Core](https://nf-co.re) is a community-driven platform that collects and shares curated sets of analysis pipelines.

#### Resources ####
* [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html)
* [Nextflow trainings ](https://nf-co.re/usage/nextflow)
* [Nextflow step by step tutorial](https://carpentries-incubator.github.io/workflows-nextflow/)

## How do I get set up? ##

#### Code Editor Installation ####
* Recommend to use [atom](https://atom.io) text editor
  * Download atom and install locally
  * Install from `Preferences > Packages` the `language-nextflow` package to enable nextflow syntax highlighting

#### Nextflow Installation ####

* [Install java](https://www.oracle.com/java/technologies/downloads/)
  * Check version `java --version` (Java 17.x.x)
    * If you don't get the right version, check its in correct `sys dir`
      * Type `which java`. On OSX this should be in `/usr/bin/java`

  * On OSX, add the below to your `.zhrc`  

  ```bash
# Check local/bin dependencies before system system/bin dependencies
PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$PATH"
  ```

* [Install nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
  * `curl -s https://get.nextflow.io | bash`
  * `chmod +x nextflow` to make `nextflow` a executable file  
  * Move the `nextflow` executable to your `$PATH`, on OSX, `mv nextflow ~/bin`


* Check `nextflow` installation
  * `nextflow -info `, this should produce a similar output to:

  ```
  Version: 21.10.6 build 5660
  Created: 21-12-2021 16:55 UTC (11:55 EDT)
  System: Mac OS X 12.1
  Runtime: Groovy 3.0.9 on Java HotSpot(TM) 64-Bit Server VM 17.0.1+12-LTS-39
  Encoding: UTF-8 (UTF-8)
  ```

## Running nf-NeST pipeline ##
Head over to [nf-NeST](https://github.com/CDCgov/Nf-NeST) and try to run it.

[To main page](https://github.com/CDCgov/MaRS) :back:  

> [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)
