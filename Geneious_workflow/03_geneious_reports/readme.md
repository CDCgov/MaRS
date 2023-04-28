> Author: @ DP updated 4/27/23 :Tiger:  
> version 0.1

----
>#### TODO ####
>#### Activity Name ####
- [ ] Add readme to each of the steps @ET :goat:
- [ ] Add gitter link @ET :goat:  

>#### Completed Activity ✓ ####
>- [x] Add Readme file: 
------

# Summarizing and Visualizing the Geneious Results #

## Background ##

This step uses the [Geneious_output_summary_and_visualization.ipynb](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_geneious_reports/Geneious_output_summary_and_visualization.ipynb) file which is writen in a python launguage, to process the Geneious output file "Annotations.csv" from step 02 [Geneious_analysis](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/02_geneious_analysis). It creates a simplyfy reports for Paired and Pooled samples separatly. It also performs Quality control and data visualization steps for the Paired and Pooled samples. 


## Requirements ##

* User first need to install the four python packages (pandas, numpy, seaborn and matplotlib) using conda to run the Jupyter code
   - conda install pandas
   - conda install numpy
   - conda install seaborn
   - conda install matplotlib
   
* [VariantsOfInterest.csv](https://github.com/CDCgov/MaRS/tree/master/Geneious_workflow/03_geneious_reports/VariantsOfInterest.csv) file which lists all the SNPs intersted that are known with drug resistance or previously reported 

* Annotations.csv file from Step02 Geneious analysis

**Note**: Save the Annotations.csv and VariantsOfInterest.csv file in same folder as the Jupyter code Geneious_output_summary_and_visualization.ipynb


## Procedure ##

* The first step users need to edit the code based on their dataset.First edit the name of csv file if it not named as “Annotations.csv”. 
* Next, the user needs to edit the site information based on the state or province of the corresponding data.  For example, as given in the Jupyter notebook, the two letter code in AMDID from digit 4:6 is AM, and it represents site name ‘Amhara’. If you have more sites, copy and paste the if statement highlighted below and add more rows depending on the number of sites present in the data being annotated.

<span style="color:red">
     def site(row):                               
         if row['Document Name'][4:6]=="Am": 
           return 'Amhara'

</span> 

* Edit the name of the output files to make it specific to the country that the data corresponds to. For example, if the country is Burkina Faso, then default output file should be changed from “ET_individual_EPI_VOI.csv” to “BF_individual_EPI_VOI.csv”. Rename all other output files to their specific country code digits.  

* Once all the modifications are made, the user can proceed to running the Jupyter code by pressing the RUN/PLAY button. It may take a few second or a few minutes for the code to run and the outputs to be generated. This will depend on the size of the dataset. If there are any errors while running the code, the user will get an error message below the code. There are some errors that may appear but will not impact the generation of the output files. It is easy to discern this by opening the output files and checking if the file is populated with information or empty.

## Output ##

* Once the code finishes running, different outputs will be generated (* denotes a wild card that will be a two letter code specific to the country being analyzed):
     * *_individual_EPI.csv: reports all the SNPs that were called by the Geneious workflow for individually sequenced samples.
     
     * *_individual_EPI_VOI.csv: reports all the SNPs of interest, as defined by the VariantsOfInterest.csv file, that were called by the Geneious workflow for individually sequenced samples.
     * *_POOLED_VOI.csv: reports all the SNPs that were called by the Geneious workflow for samples that were sequenced as pools.
     * *_POOLED_EPI_VOI.csv: reports all the SNPs of interest, as defined by the VariantsOfInterest.csv file, that were called by the Geneious workflow for samples that were sequenced as pools.
     * *_weighted_bysite_EPI.csv: reports a weighted variants allele frequency(VAF) average for each SNPs by site/province
     * *_Tab_Table_SNP_Total1_Comb.csv: reports each gene with its Major (AF >= 50%), Minor (AF < 50%) and Wildtype (AF = 0%) allele frequency. 
     * *_Bar_plot_Combined3: Bar graph shows the wild type, major and minor allele frequencies of all the reportable SNPs described in VariantsOfInterest.csv file. Allele frequencies are indicated on the x axis, and the variants of interest and total number of samples are listed along the y-axis.          The color coding indicates the type of mutation found in the samples; blue is for wild type, green for minor allele mutation and red for major allele mutation.

