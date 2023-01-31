> [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### Background ###

The AMD genomics lab uses a standardized samples naming schema, called `AMD ID`, to capture sample associated meta-data prior to sample prep and sequencing. By creating the `AMD ID` for each sample at this stage, the chances of erroneous linking of information post-sequencing are decreased and need to use multiple different databases for linking information eliminated. We highly encourage to use the [AMD create template](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/files/AMD_ID_create_template.xlsx) to generate the `AMD ID` prior to starting any experiments. For additional details regarding the `AMD ID` see this [presentation](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/files/AMD_ID_create_key.pptx).

In this part of the workflow, a `jupyter-lab notebook` will be used to check the `AMD IDs` of the received `AMD_ID.fastq` files. Before we start, however, let's review the two sequencing methods used by the AMD Genomics lab.  

The lab currently uses two sequencing methods: `individual` and `pooled`.

`Individual` samples are tagged using unique indices prior to sequencing to allow de-multiplexing at the individual sample level. The individual samples often have more than one molecular marker included (e.g. pooled molecular markers).

The `pooled` samples are combined (e.g. pooling of individual samples) either as filter blot spots or DNA  based on parasitemia prior to molecular marker amplification and then pooling (e.g. pooling of molecular markers per sample pool). See recent publication for more info [here](https://pubmed.ncbi.nlm.nih.gov/35030215/). Below is an example `AMD ID` for each sequencing method:

<img
  src="/images/ind_sample_ID.png"
  width="600">

  <img
    src="/images/pooled_sample_ID.png"
    width="600">


### How do I get set up? ###
* If you haven't already done so, set up your working environment for `jupyater-lab`:
  * [Virtual and python environment](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/py_ve_setup.md)
  * Additional required dependencies and modules will be listed in each jupyter-lab notebook
    * These you install using the `jupyter-lab` terminal by `python3 -m pip install <module name>`


* Clone a local copy of the **MaRS project**: `git clone https://github.com/CDCgov/MaRS.git`
* Example data can be downloaded using the [data_download jupyter notebook](https://github.com/CDCgov/MaRS/blob/goat_dev/Geneious_workflow/Test_data/data_download.ipynb) in the [Test data directory](https://github.com/CDCgov/MaRS/tree/goat_dev/Geneious_workflow/Test_data); please follow directions provided in the jupyter-lab notebook, including the steps to convert the SRAID to AMDID.  
* To run a test, start a `jupyter-lab notebook` in your terminal `cd` to location of downloaded **MaRS project directory**
  * Then run: `jupyter-lab MaRS/`
    * This will open up `jupyter-lab notebok` and show the full directory structure of the **MaRS project**
* Next navigate to: `Geneious_workflow/01_sample_ID_QC` and open `AMDID_QC_fastq.ipynb`


----
<details>
  <summary><strong>TODO</strong></summary>


> Author: :goat: ET @ 04/5/22
>> Edited & Reviewed::goat: @ET 11/30/22 
----
>#### TODO ####
>#### Activity Name ####

>#### Completed Activity ✓ ####
  - [x] Update readme to point to new example fastq (external) @ET
   - [x] Add directions on how to download from SRA @DP

- [x] Update readme @ET
    - [x] add info on naming schema for AMD_IDs  
    - [x] instructions on how to generate AMD_IDs
    - [x] test links to images  
 - [x] Jupyter-lab set up @ET 
    - [x] How to install?
    - [x] How to start up, run quick test
 - [x] Add example input data

------
</details>
