# README #

[Learn README formatting using Markdown](https://bitbucket.org/tutorials/markdowndemo)

## What is this repo for? ##

* This README will document all the steps necessary to get [nextflow](https://www.nextflow.io/docs/latest/index.html)
up and running using the [NeST](https://github.com/CDCgov/NeST) pipeline, including how to maintain and further
develop it.
* Version: 0.1

#### Nextflow ####

[Nextflow](https://www.nextflow.io/docs/latest/index.html) is a bioinformatics workflow management software that enables scalable and reproducible workflows. It works with Docker, Singularity and Conda. It uses a common scripting language and processes, allowing for easy integration and running on cloud or on-premise HPC environments. [Nf-Core](https://nf-co.re) is a community-driven platform that collects and shares curated sets of analysis pipelines.

#### Resources ####
* [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html)
* [Nextflow trainings ](https://nf-co.re/usage/nextflow)

## How do I get set up? ##
* [Installation of code editor](#Code_Editor)
* [Installation of nextflow ](#Nextflow)

<a id="Code_Editor"></a>
## Code Editor Installation ##
* Recommend to use [atom](https://atom.io) text editor
  * Download atom and install locally
  * Install from `Preferences > Packages` the `language-nextflow` package to enable nextflow syntax highlighting

<a id="Nextflow"></a>
## Nextflow Installation ##

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

<a id="Nextflow Script"></a>
## Run an example Nextflow script ##
Below is a well documented nextflow script that takes as input `fastq.gz` files and counts the number
of lines in each `fastq.gz`

```java
{nextflow codeblock}
```

<a id="Nextflow Processes"></a>
## Nextflow Tutorial ##

See [Nextflow step by step tutorial](https://carpentries-incubator.github.io/workflows-nextflow/).

To control inputs, outputs and how a command is executed a process may contain five definition blocks:

- `directives`: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses and the amount of memory allocated.
- `inputs`: Define the input dependencies, usually channels, which determines the number of times a process is executed.
- `outputs`: Defines the output channels used by the process to send results/data produced by the process.
- `when clause`: Allows you to define a condition that must be verified in order to execute the process.
- `script block`: A statement within quotes that defines the commands that are executed by the process to carry out its task.

The syntax is defined as follows:
```java
    process < NAME > {
     [ directives ]        
     input:                
     < process inputs >
     output:               
     < process outputs >
     when:                 
     < condition >
     [script|shell|exec]:  
     < user script to be executed >
    }
```
The `input` qualifier declares the type of data to be received. Types of  qualifiers:

- `val`: Lets you access the received input value by its name in the process script.
- `env`: Lets you use the received value to set an environment variable named as > the specified input name.
- `path`: Lets you handle the received value as a path, staging the file properly in the execution context.
- `stdin`: Lets you forward the received value to the process stdin special file.
- `tuple`: Lets you handle a group of input values having one of the above qualifiers.
each: Lets you execute the process for each entry in the input collection.





* [Placeholder](#Placeholder)
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
