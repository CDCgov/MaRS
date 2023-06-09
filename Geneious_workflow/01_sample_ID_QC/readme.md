1qwr4> Author: @ET 4/5/22 :goat:  
>> Edited: @DP 4/25/22: Tiger
----
>#### TODO ####
>#### Activity Name ####

>#### Completed Activity ✓ ####
 - [x] Update readme
    - [x] add info on naming schema for AMD_IDs  
    - [x] instructions on how to generate AMD_IDs
    - [x] test links to images  
 - [x] Jupyter-lab set up
    - [x] How to install?
    - [x] How to start up, run quick test
 - [x] Add example input data

------

### Background ###

The AMD (Advanced Molecular Detection) genomics lab uses a standardized samples naming schema, called `AMD ID`, to capture sample associated meta-data prior to sample prep and sequencing. By creating the `AMD ID` for each sample at this stage, the chances of erroneous linking of information post-sequencing are decreased and need to use multiple different databases for linking information eliminated. We highly encourage to use the [AMD create template](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/files/AMD_ID_create_template.xlsx) to generate the `AMD ID` prior to starting any experiments. For additional details regarding the `AMD ID` see this [presentation](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/01_sample_ID_QC/files/AMD_ID_create_key.pptx).

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

#### AMD ID Description ####

* The AMD ID Key:  `<year> <country> <state/site> <day of failure> <treatment> <sample_id> <genus spp> <sample type> <mol marker bit code> <# sample processed>`.
* Any missing meta data is replaced with an `x` _lower case_ strings for each character position.

Example:
- `Individual` sequenced sample ID: `17GNDo00F0001PfF1290` = `<2017> <Guinea> <Dorota> <Day0> <AS+AQ> <0001> <P.falciparum> <FilterBloodSpot> <k13-crt-mdr-dhfr-dhps-cytB-cpmp-pfs47> <0>`

- `Pooled` sequenced sample ID: `17GNDo00x001p10F1290` = `<2017> <Guinea> <Dorota> <Day0> <missing info> <001> <Pooled Samples> <Samples in Pool> <FilterBloodSpot> <k13-crt-mdr-dhfr-dhps-cytB-cpmp-pfs47> <0>`

NOTE: If information is not availble (na), **x** is used for each character position. For example, in the pooled samples Treatment has (2) character spaces, represented as a two digit number. This is replaced with (2) **xx** since its a pool of samples that have possibly different treatment information. Moreover, for pooled samples, **Sample ID** is replaced with **three digit number + the letter p** (for pooled), and **Genus** is replaced with **total number of SamplesInPool** as a (2) digit number.



### How do I get set up? ###
* If you haven't already done so, set up your working environment for `jupyater-lab`:
  * [Virtual and python environment](https://github.com/CDCgov/MaRS/blob/master/Geneious_workflow/py_ve_setup.md)
  * Additional required dependencies and modules will be listed in each jupyter-lab notebook
    * These you install using the `jupyter-lab` terminal by `python3 -m pip install <module name>`


* Clone a local copy of the **MaRS project**: `git clone https://github.com/CDCgov/MaRS.git`
* Test data will need to be downloaded using the data_download jupyter notebook in the Test data directory
* To run a test, start a `jupyter-lab notebook` in your terminal `cd` to location of downloaded **MaRS project directory**
  * Then run: `jupyter-lab MaRS/`
    * This will open up `jupyter-lab notebok` and show the full directory structure of the **MaRS project**
* Next navigate to: `Geneious_workflow/01_sample_ID_QC` and open `AMDID_QC.ipynb`

### Procedure  ###

* Using jupyter-lab notebook (AMDID_QC.ipynb), the investigator should check the validity of the newly generated AMD IDs found in named fastq files. This Jupyter code will check each sample ID for each fatsq file or sample ID listed in csv file(Samplesheet) before samples proceed to sequances to determine if it adheres to the standardized AMD ID format. It will confirm the character count and ensure each ID contains the required elements described in AMD ID description.

* Once jupyter-lab is launched, it will open in the web browser (Chrome, Edge, Firefox etc.). Before running the quality check Jupyter Notebook, the user also needs to install specific python packages (Pandas, tabulate) to run the jupyter code. This can be done by opening terminal from inside the Jupyter-lab and install the packages using conda as shown below.
  * conda install pandas
  * conda install tabulate

* The quality chek code and documentation is contained within the file called AMDID_QC.ipynb. Open this file and run the each code section by pressing the play button or pressing the control/command and Enter.

* User first needs to edit the code by providing the full path of the folder where the fastq files are located (**/User/admin/Desktop/Folder_with_fastq/.fastq.gz**)
* If User needs to veryfy the sampleIDs from a samplesheet even before runing the sequencer, provide the csv path of that samplesheet (**/User/admin/Desktop/sample.csv**)

* After running the code, if there are no errors in the sample naming, the user will get the message to proceed to next step as all the samples passe through the QC test. If there is an error in sample naming, the user will get the list of sample IDs which did not match with the standardized AMD IDs. Ensure the user makes the necessary corrections before proceeding.
