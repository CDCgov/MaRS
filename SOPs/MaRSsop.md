# Targeted Amplicon Deep Sequencing of *Plasmodium falciparum (Pf)* molecular markers of resistance

*Preparing full length genes Pf kelch 13 (k13), Pf chloroquine resistant transporter (crt), Pf multidrug resistant protein 1 (mdr1), Pf bifunctional dihydrofolate reductase-thymidylate synthase (dhfr), Pf dihydropteroate synthase (dhps) and mitochondrial genome amplicons for sequencing on the Illumina MiSeq system*

*If you end up using MaRS in your workflow, please cite this [study](https://www.ncbi.nlm.nih.gov/pubmed/29439965):*

```
Next-Generation Sequencing and Bioinformatics Protocol for Malaria Drug Resistance Marker Surveillance.

Talundzic E, Ravishankar S, Kelley J, Patel D, Plucinski M, Schmedes S, Ljolje D, Clemons B, 
Madison-Antenucci S, Arguin PM, Lucchi NW, Vannberg F, Udhayakumar V.

Antimicrob Agents Chemother. 2018 Mar 27;62(4). pii: e02474-17. doi: 10.1128/AAC.02474-17. Print 2018 Apr.
```
---------

# Table of Contents
 * [Introduction/General Overview](#chapter-1)
 * [Materials and Equipment](#chapter-2)
 * [NGS Library Preparation Workflow](#chapter-3)
 * [PET PCR Sample QC](#chapter-4)
 * [Gene PCR Enrichment](#chapter-5)
 * [Electrophoresis](#chapter-6)
 * [SequalPrep Normalization and Purification](#chapter-7)
 * [Tagment Genomic DNA](#chapter-8)
 * [Library Amplification and Index PCR](#chapter-9)
 * [Library Clean up](#chapter-10)
 * [Next Steps](#chapter-11)
 * [Library Pooling, Quantification, and Normalization](#chapter-12)
 * [Library Denaturing and MiSeq Sample Loading](#chapter-13)
 * [Supporting Information](#chapter-14)

-----------

**IMPORTANT NOTICE** 
This document provides information for an application for Illumina technology that has been demonstrated internally and may be of interest to customers. This information is provided as‐is and is not an Illumina product and is not accompanied by any rights or warranties. Customers using or adapting this information should obtain any licenses required and materials from authorized vendors. Illumina products mentioned herein are for research use only unless marked otherwise. While customer feedback is welcomed, this application is not supported by Illumina Technical Support and Field Application Scientists.

November 2nd, 2017 version 3.3                          modified by Eldin Talundzic (etalundzic@cdc.gov)
This protocol format was adapted from the 16S Metagenomics Sequencing Library Preparation protocol from Illumina.

-----------

<a id="chapter-1"></a>
## Introduction and General Overview 
Human malaria is caused by four parasite species called *Plasmodium falciparum (Pf), P. vivax (Pv), P. ovale (Po)* and *P. malariae (Pm)*. Two of these, Pf and Pv, pose the greatest threat to global public health. About 3.2 billion people, half of the world’s population, are at risk for malaria. In 2015, 95 nations had ongoing malaria transmission and annually 198 million people suffer from malaria, causing an estimated 438,000 deaths. In the U.S., an estimated 1,500 - 2,000 cases of malaria are imported annually. One of the greatest public health challenges for malaria control and elimination is the threat of drug resistant Pf parasites. 

Our previously effective anti-malarial treatments, chloroquine (CQ) and sulfadoxine/pyrimethamine (SP), are ineffective in many regions. Even more alarming, resistance to our last effective class of anti-malarial drugs, called artemisinins, has now emerged and spread in Southeast Asia — threatening malaria control and prevention programs globally. 

Identifying and tracking drug resistance, are critical for providing appropriate malaria prophylaxis and treatment guidelines. Molecular markers of resistance are available for several anti-malarial drugs, including artemisinins. Surveillance using molecular markers provides a robust system for the detection and tracking of resistant malaria parasites.
Below is a table of the major molecular markers and their respective amino acid polymorphisms that are associated with resistance to various anti-malarial treatments. 


**Table 1** Molecular Markers and amino acid polymorphisms

**_Resistance to treatment_** | **_Molecular marker_** | **_Associated drug resistance mutations_**
|------------------------ | --------------------------- | -----------------------------------------
*Chloroquine* | *crt* | M74I    N75E    K76T    A220S    Q271E    N326S    I356T    **\*C350R**    R371I 
*Chloroquine, quinolone, mefloquine, aminoquinoline halofantrine* | *mdr1* | N86Y    Y186F    S1034C    N1042D    D1246Y 
*Atvaqoune* | *cytochrome b* | I258M    Y268S
*Artesunate Sulfadoxine-pyrimethamine* | *kelch 13* | N458Y    Y493H    R539T    I543T    R561H    C580Y 
*Artesunate Sulfadoxine-pyrimethamine* | *dhps* | S436A/F    G437A    K540E/G540E    A581G    A613S/T
*Sulfadoxine-pyrimethamine* | *dhfr* | N51I    C59R    S108N

_**\*C350R** results in a reduction of chloroquine resistance._

The method described in this SOP provides an integrated advanced molecular detection system for drug resistance surveillance of all currently known anti-malarial treatments. This technology is based on a next generation sequencing (NGS) protocol referred to as targeted amplicon deep sequencing (TADS). 

<a id="chapter-2"></a>
## Materials and Equipment 
Please ensure all the necessary user‐supplied consumables and equipment are available before proceeding to sample preparation.


**Table 2** User‐Supplied Consumables

Consumable | Supplier
------------------------------------- | --------------------------------------------
**Non-powdered sterile gloves** | General lab supplier
**Laboratory coat** | General lab supplier
**1.7ml microcentrifuge tubes** | General lab supplier
**10 μl barrier pipette tips** | General lab supplier
**10 μl multichannel pipettes** | General lab supplier
**10 μl single channel pipettes** | General lab supplier
**20 μl barrier pipette tips** | General lab supplier
**20 μl multichannel pipettes** | General lab supplier
**20 μl single channel pipettes** | General lab supplier
**200 μl barrier pipette tips** | General lab supplier
**200 μl multichannel pipettes** | General lab supplier
**200 μl single channel pipettes** | General lab supplier
**1000 μl barrier pipette tips** | General lab supplier
**1000 μl multichannel pipettes** | General lab supplier
**1000 μl single channel pipettes** | General lab supplier
**PCR grade water** | General lab supplier
**RNase/DNase‐free 8‐well PCR strip tubes and caps** | General lab supplier
**\[Optional] Disposable Polystyrene Reservoirs** | General lab supplier (Thomas Scientific Cat # 55501080)
**2X ABI TaqMan environmental buffer w/ Rox dye** | Applied Biosystems Catalog # 4396838
**Strip tubes 8X** | Agilent Catalog #410022
**Strip Optical caps 8X** | Agilent Catalog #410024
**Hardshell® 96‐well PCR Plates, clear semi skirted** | Bio-Rad Catalog # HSS9601 
**Phusion® High-Fidelity DNA Polymerase** | NEB Catalog # M0530L
**Deoxynucleotide (dNTP) Solution Mix** | NEB Catalog #N0447L
**Amplicon PCR Forward Primer (Standard desalting)**	| (See tables 2-7)
**Amplicon PCR Reverse Primer (Standard desalting)** | (See tables 2-7)
**Lonza SeaKem® LE Agarose** | Lonza Catalog # 50004
**Nucleic Acid gel stain** | Biotum GelRed™ Nucleic Acid Gel Stain
**DNA gel loading dye** | Yakva Scientific 6X Orange-G Gel Loading Buffer  #YSG Lot #36
**Quick-Load 1kb DNA ladder** | NEB Catalog # N0468
**UltraPure™ 10X TBE Buffer** | Fisher Scientific Catalog #15581-044
**SequalPrepTMNormalization Plate Kit, 96-well** | ThermoFisher Scientific Catalog # A1051001
**Nextera® XT DNA Sample Preparation Kit (96 Samples)** | Illumina, Cat # FC-131-1096
__\* Nextera XT Index Kit (96 Indexes, 384 samples)__ | Illumina, Cat # FC‐131‐1002
**Agencourt AMPure XP 60 ml kit** | Beckmen Coulter Genomics, part # A63881
**200 Proof Ethanol** | Decon Labs, Inc. Catalog # 2716
**AlumaSeal II aluminum seals** | Excel Scientific, Inc. Cat # AF100
**Clear, 8-strip PCR tubes domed caps** | LabSource, Item # T54-203-CS/10PK MFG# - 321-10-062
**\[Optional] 96‐well storage plates, round well, 0.8 ml (“MIDI” plate)** | Fisher Scientific, part # AB‐0859
**Qubit® dsDNA HS Assay Kit** | Life Technologies Corporation Catalog #Q32853
**Qubit™ Assay Tubes** | Thermo Fisher Scientific Catalog # Q32856
**Agilent High Sensitivity D5000 ScreenTape** | Agilent Technologies, Cat. # 5067-4626 (5067-5592)
**Agilent High Sensitivity D5000 Reagents** | Agilent Technologies, Cat. # 5067-4626 (5067-5593)  
**Agilent High Sensitivity D5000 Ladder** | Agilent Technologies, Cat. # 5067-4626 (5067-5594) 
**Quant-iT™ PicoGreen® dsDNA Assay Kit** | Life Technologies, Cat # P11496
**1N NaOH** | Sigma-Aldrich, Inc. Catalog #SX0607H-6
**Tris Hydrochloride, 1M Solution (pH 7.0/Mol. Biol.)** | Thermo Fisher Scientific Catalog #BP1756-100 
**MiSeq Reagent Kit V2 500 cycle kit** | Illumina #MS-102-2003
**MiSeq Reagent V2 Nano Kit 500 cycle** | Illumina #MS-103-1003

-----------------

__* If you plan to pool >96 samples, you will need the Nextera XT Index Kit v2 Set A, B, C, and D to provide unique multiplex combinations of indices.__


**Table 3** User‐Supplied Equipment

Equipment | Supplier
-------------------------------- | ------------------------------------
**2-8°C Refrigerator** | General lab supplier
**-20°C Refrigerator** | General lab supplier
**Vortex** | General lab supplier
**4x Eppendorf PCR Cooler, iceless cold storage system for 96 well plates and PCR tubes** | Sigma-Aldrich, Inc (Z606634-1EA)
**Agilent ABI7500 or equivalent real-time PCR macine** | Agilent Technologies, Cat # 4351106
**96‐well thermal cycler (with heated lid)** | General lab supplier
**Electrophoresis rig** | General lab supplier
**Magnetic stand‐96 or 96S Super Magnet Plate** | Life Technologies, Cat # AM10027 or Alpaqua SKU A001322 
**Microplate centrifuge** | General lab supplier 
**Qubit 3.0 Fluorometer** | Life Technologies Corporation, Cat #Q33216
**Agilent D4200 ScreenTape System** | Agilent Technologies, Cat # G2991AA
**MiSeq Desktop Sequencer** | Illumina Inc.

<a id="chapter-3"></a>
## Protocol Workflow 
**NOTE: The hands-on times are based on using 96-well format plates for each step.**

#### PET-PCR Sample Quality Check \[Sample QC]
Real-time PCR hands-on time 30 min / 96 samples; Cycle time 1.2 hours  
Reagents: Primers, 2X ABI TaqMan buffer, DNase PCR free water  

#### PCR reaction to generate amplicons \[Amplification]
PCR hands-on time 30 min / 96 samples; Cycle time 2.5 hours  
Reagents: 10uM Primers, HF Phusion Taq, 5X GC Buffer, 10mM dNTPs, DNase PCR free water  

#### Analysis of PCR amplicons \[Electrophoresis]
PCR amplicon electrophoresis hands-on time 10 min / 8 samples; Gel running time 30 min  
If < 20 samples, run all samples on the gel; If > 20 samples, pick 20 samples with varying CT values and run on the gel   
Reagents: Agarose, DNA loading dye, 1kb DNA ladder, 1X TBE Buffer  

#### PCR amplicons clean up \[Purification]
Hands on time 35 min / 96 samples; Total time 90 min / 8 samples  
Reagents: SequalPrep Normalization Binding Buffer, SequalPrep Normalization Wash, SequalPrep Normalization Elution Buffer  

#### Tagment Genomic DNA
Hands on time 30 min / 96 samples; Total time 17 min / 8 samples  
Reagents: ATM, TD, NT  

*\[optional]* To assess tagmentation, run 1μL sample on Agilent Bioanalyzer 2X and/or TapeStation 2X using High Sensitivity DNA chip  

#### Library Amplification
Hands on time 35 min / 96 samples; Cycle time 38 min / 96 samples  
Reagents: NPM, Index 1 primers, Index 2 primers  

#### Library Clean-up \[Purification]
Hands on time 30 min / 96 samples; Total time 40+ min / 96 samples  
Reagents: AMPure XP beads, fresh 80% EtOH  

#### Library Pooling, Quantification, and Normalization
Hands on time 30+ min / 96 samples; Total time 40+ min / 96 samples  
Reagents: Sample Buffer, D5000 Ladder, ScreenTape; Qubit dsDNA HS Buffer and Reagent, Standard #1 and #2  

#### Library Denaturing and MiSeq Sample Loading
Hands on time 30 min / pooled samples; Total time 30 min / pooled samples  
Reagents: Resuspension Buffer, HT1, 0.2N NaOH, PhiX Control Kit v3, 200mM Tris-HCl pH7.0  

#### Analysis of NGS data \[Analysis]
Hands on time 5 min / 96 samples; Total time 15-25 min / 96 samples  
Method: MaRS analysis pipeline  
__*Standardized SNPs reports generated*__  

<a id="chapter-4"></a>
## PET-PCR Sample QC 
This step uses a real time PCR assay to assess the quality and quantity of starting DNA material. The readout includes an estimation of *all* DNA in the sample, host and parasite. 

**NOTE:** This assay is used to identify which samples should be included and/or excluded for downstream procedures. Any sample with a CT value greater than 34, should be excluded. We recommend grouping all samples with a CT > 34 and performing separate amplicon PCRs and electrophoresis for these samples. If the electrophoresis yields positive results (e.g. positive and clear bands on gel) you may procced with downstream procedures. Samples with a CT > 34 have either very low and/or no parasite DNA.

### Consumables

**Table 4. PET-PCR Consumables**

Item | Quantity | Storage
------------------------------- | ------------------ | ----------------
**Primers – FAM labled genus primers and HEX labeled *falciparum* primers (see below)** | 0.25-0.5 μl per sample | 2° to 8°C
**TaqMan 2X Environmental buffer** | 10 μl per sample | 2° to 8°C
**Nuclease-free water** | 6.25 μl per sample | Room temperature
**Strip tubes 8X** | Up to 8 samples per strip | Room temperature
**Strip Optical caps** | 8X	Up to 8 samples per strip | Room temperature

### Preparation

* All stock primers should be prepared at a 10μM concentration.
* DNA samples should be stored at 4°C until testing or -20°C for long term storage.
* Store all primer stocks at -20°C for up to 1 year.
* Unopened tubes of ABI TaqMan Environmental Buffer should be stored at -20°C for a maximum of six months. Once thawed, store at 4°C for up to six months. The reagent must be used within the expiration date provided by the manufacturer.
* All samples should be tested in duplicates or triplicates in some special cases (e.g. very low density situations).
* Typically, the genus/*P. falciparum* multiplex assay should be run first on all samples. All genus-positive sample are subsequently tested for *P. ovale*, *P. malaria* and *P. vivax*, in order to determine the species.
    
### Procedure

**Initial Set up**

* The PET-PCR reaction mix is prepared by mixing the TaqMan environmental buffer, primers, and water as shown below. 
* Determine the number of reactions you need to run by multiplying the total number of samples you have to test (including your positive and negative controls) by two because every sample will be tested in duplicates.
    * For example, if you are testing 10 samples, you will multiple this by two to give you 20. Add two extra reactions to account for loss of solution during pipetting. This gives you a total number of 22 reactions.  Multiply this number with the volumes below for each component to give you the total master-mix volume required for your experiment. 
* In a 1.5mL tube, prepare your master-mix by multiplying the volumes shown below with the total number of reactions you need to run (e.g. 22 as described above). 

#### Primers and PCR Conditions

The table below shows the primers and PCR conditions for a multiplex reaction- Genus and *P. falciparum*:

**Table 5. _Multiplexing genus and P. falciparum species specific primers_** 

|                               | _**20.0 μl rxn**_ | X samples +1 | _**Final \[conc]**_ |
|-------------------------------- | ----------------- | ------------- | -------------------|
| *H<sub>2</sub>O* | 6.25 μl | | |
| *2X ABI TaqMan buffer* | 10.0 μl | | *1 x*|
| *Genus F primer (10uM)* |0.5 μl |  | *250 nM*|
| *FAM-genus R primer (10uM)* | 0.5 μl |  | *250 nM*|
| *P. falciparum F primer (10uM)* | 0.5 μl |    | *250 nM*|
| *HEX-P falciparum R primer (10uM)* | 0.25 μl |    | *125 nM*|
| *Total* |**18.0 μl** |   |    |                         
| __*Add last:*__ |     |    ||
| *DNA* |2.0 μl |   ||
| **TOTAL** |20.0 μl | ||

--------------------------
   
|      | *Thermal cycling conditions* |  |
| ----: | :------------: | --- |
|*95°C* |*15:00 min* | |
|__*95°C*__ | __*0:20 min*__ | |
|__*63°C*__ |__*0:40 min*__ | __*45X*__ |
|__*72°C*__ | __*0:30 min*__ | |

**Primers:**  
**Genus 18sFor (5' to 3'):**  5’-GGC CTA ACA TGG CTA TGA CG-3’
**Genus 18sRev (5' to 3'):** 5’-agg cgc ata gcg cct ggC TGC CTT CCT TAG ATG TGG TAG CT-3’ **(FAM-labeled: based on the 18s rRNA gene)**

**P. falciparum For (5' to 3'):** 5’-ACC CCT CGC CTG GTG TTT TT-3’ 
**P. falciparum Rev (5' to 3'):** 5’-agg cgg ata ccg cct ggT CGG GCC CCA AAA ATA GGA A-3’ **(HEX-labeled: based on the r364 target)**


#### Adding the DNA Samples

1. Mix the prepared master-mix well by vortexing briefly.
2. Centrifuge the tubes for 5 seconds to remove any solution trapped in the cap.
3. Arrange the optically clear PCR tubes on a PCR-tube rack following the PCR sample sheet. Add 18 μl of the PET-PCR master mix prepared above to each PCR well. Loosely put on the lids of the wells filled with master mix solution.
4. Return all reagents to the freezer and refrigerator before proceeding to the next step.
5. Take the assembled plate containing the tubes with PCR master mix solution to the PCR template area. 
6. Add 2μl of the unknown DNA samples to the wells with the master-mix according to the sample sheet. Cap the well tightly after adding the sample. The total volume of PCR reaction is 20.0 μl after addition of the template.
7. Add positive control DNA to each positive control well with master-mix. Cap the wells after each positive control is added. 
8. Add 2.0 μl of DNase-free H<sub>2</sub>O to the wells designated as the no-template control (NTC) and close that well tightly.
9. Make sure each sample has been added to the correct well and that all wells are tightly capped.
10. Briefly centrifuge your strip tubes to remove any solution trapped on the walls of the wells.
11. Make sure there are no bubbles in the well.

**NOTE:** *The amount of template DNA to be used can be as low as 2μL, but it is not uncommon to use 5μL.  This can be adjusted appropriately depending on the sample parasitemia. The change should be discussed before it is implemented.*

#### PCR-Cycling Parameters

1. Start the real-time PCR thermocycler according to the manufacturer’s guidelines.
2. Program the software to detect fluorescence through FAM, HEX and ROX filters all wells. ROX is to be detected as a reference dye.
3. Program the software to run the cycling conditions shown on page 7.
4. Fluorescence data should be collected at the amplification plateau.

#### Interpreting Results

1. Interpret the results using standard settings in the software.
2. If the calculated thresholds are located within the background noise, they should be manually set to a level slightly higher than the background. Such alterations should be done with only one dye displayed at the time.
3. Positive specimens are those that yield a fluorescence signal above the threshold value in the wells where samples or controls were loaded
    * Positive PCR: A positive sample produces a fluorescence signal above the threshold/noise level. Positive samples are designated a Ct value below 40.0.
    * Negative PCR: No fluorescence signal above the threshold/noise level. Negative samples have no Ct or have a Ct value above 40.0.
 
**NOTE:** *The negative controls must be negative (no Ct or above 40.0). The positive controls must be positive (designated by Ct value below 40.0). The test should be repeated if the NTC has a positive Ct value, or if the positive control yields no positive results.*

**For more information, please see:** 
Lucchi, N.W., et al., Molecular diagnosis of malaria by photo-induced electron transfer fluorogenic primers: PET-PCR. *PLoS One*, 2013. 8(2): p. e56677.

<a id="chapter-5"></a>
## Gene PCR Enrichment 
##### This step uses PCR to amplify template from a DNA sample using region of interest-specific primers.
User‐defined forward and reverse primers are used to amplify templates from genomic DNA. A subsequent limited‐cycle amplification step is performed to add multiplexing indices and Illumina sequencing adapters. Libraries are normalized and pooled, and sequenced on the MiSeq system using v2 reagents.

#### Procedure 
**Initial Set up**
- Ensure that the No-DNA and DNA-only UV stations have all the appropriate pipettes and tip sizes.
- Clean up all pipettes and lab bench area using using 10% bleach followed by 70% ethanol.
- Turn on the No-DNA (PCR master mix) and DNA-only UV station for 30 minutes.  
- If you have not already done so, create a 10uM working stock solution of your primers, using C1V1 = C2V2 to calculate the appropriate volume needed to make a working stock solution. 
- Label all freshly made and newly opened items with the date and your name initials.
- Get the appropriate number of PCR plates and/or PCR tubes and place them in the no-DNA UV station. 
- Label the PCR tubes (use a printout template for PCR plates) with sample IDs. 
- Get the appropriate number Eppendorf PCR Cooler plates from the freezer, wipe down with 70% ethanol and place them in the UV hoods (one or more sets each in the no-DNA and DNA-only UV stations). Turn on the UV stations for another 30 min.
- Be sure to reserve the appropriate number of thermocyclers and have the appropriate cycling conditions set up (**Table 6**).

**Step by step procedure**
- Let Primers, dNTPs, and GC Buffer defrost at room temperature (10-15 min). Once defrosted, mix gently (vortexed) and centrifuge briefly prior to use. __*DO NOT thaw and/or vortex or mix the HF Taq.*__
- All PCR reactions **must be assembled on the** Eppendorf PCR Cooler plates
- **Always add the Taq last when making your master mix and DO NOT vortex and/or pipette after adding Taq.** 
- If you forget to return any of the reagents, especially the Taq, to its appropriate storage conditions (i.e., leave it out at room temperature), record the date and time of when it happened, and discard.  


1. Set up the following reaction of water, GC Buffer, dNTPs, primers, HF Taq Phusion, and DNA in the order given in **Tables 6.1 – 6.6**:
    - Calculate appropriate volumes for mastermix based on number of samples to be included in reaction; multiply each reagent volume times the total number of samples + 1 (for user pipetting errors)
    - Final volume of master mix is given in **Tables 6.1 - 6.6**.
    - _**NOTE:** If the number of samples is <5, make a mastermix for at least 6 samples to avoid pipetting volume  errors._
2. Seal plates and/or PCR tubes. 
3. Once tubes and/or plates are sealed, keep them in the Eppendorf PCR Cooler plates. **Pre-heat** the thermal cycler to 98°C prior to placing PCR plates and/or PCR tubes into the thermal cycler. Pre-heating to 98°C should take 0:30 of the 3:00 min. 


### Primers and PCR Conditions 
The tables below show primers and PCR conditions for **Pfcrt (6.1), Pfk13 (6.2), Mitochondrial genome (6.3), Pfcytb (6.3a), Pfdhps (6.4), Pfdhfr (6.5), Pfmdr1 (6.6).**

**Table 6.1 *Pfcrt***

|*Pfcrt (3109bp)*    | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                           
|:------------------:|:--------------:|:-------------:|:--------------:|
| *Master Mix:*      |                                                                                                 	                                                             
| *5X GC Buffer*    | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)* |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H<sub>2</sub>O*	            |  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl** |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |


|Thermal  |cycling  |conditions|
|:------:|:---------:|:----------:|
|*98°C*        |*3:00 min*  |           |
|***98°C***	|**_0:10 min_**|           |
|***62°C***	|**_0:30 min_**| **_30X_** |
|***65°C***	|**_5:00 min_**|       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   




**Primers:  Pfcrt_F1 Sequence (5' to 3'):** TTACATATAACAAAATGAAATTCGC; **Pfcrt_R1 Sequence (5' to 3'):** TATTGTGTAATAATTGAATCGACG


**Table 6.2 _Pfk13_**

|*Pfk13 (2120bp)*	  | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                           
|:------------------:|:--------------:|:-------------:|:--------------:|
| *Master Mix:*      |                                                                                                 	          | *5X GC Buffer*     | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)*  |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H<sub>2</sub>O*	 |  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl** |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |


|Thermal  |cycling  |conditions|
|:------|:---------|:----------|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***58°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   



**Primers:  Pfk13_F1 Sequence (5' to 3'):** CTATGACGTATGATAGGGAATCTGG; **Pfk13_R1 Sequence (5' to 3'):** CTGGGAACTAATAAAGATGGGCC




**Table6.3 Mitochondrial genome**

|*Mitch(5967bp)*	  | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                           
|:--------------------:|:--------------:|:-------------:|:-----------------:|
| *Master Mix:*      |                                                                                                 	                                               
| *5X GC Buffer*    | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)* |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H<sub>2</sub>O*	            |  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl** |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |
 
 |Thermal  |cycling  |conditions|
|:------|:---------|:----------|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***50°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   




**Primers: Mitochondrion_F1 Sequence (5’ to 3’):** AAGCTTTTGGTATCTCGTAAT; **Mitochondrion_R1 Sequence (5’ to 3’):** TATTATAATATAACTCTACAAAGTTGAAC




**NOTE:** If experiencing issues with amplifying the full-length mitochondrial genome, consider amplifying only the cyt-b gene instead for characterizing molecular markers associated with Malarone (atovaquone/proguanil) resistance. See Table 6.3a below. 




**Table6.3a *Cytochrome b***


|*Pfcytb 1 (937bp)*	 | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                           
|:------------------:|:--------------:|:-------------:|:--------------:|
| *Master Mix:*      |                                                                                                 	                                                       
| *5X GC Buffer*     | 	*10.0 μl*      |                  | *1 x*|
| *dNTPs (10mM)*	   |   *1 μl*	     |                  | *0.2 mM*|
| *F primer (10uM)*  |	*1.25 μl*	  |                  | *0.25 uM*|
| *R primer (10uM)*	 |  *1.25 μl*  |                  |*0.25 uM*|
| *H<sub>2</sub>O*    |  *32 μl*   |                  |         |
 |**Add last:** 	   |             |                  |         |
 |*HF Phusion Taq*   |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	           | **46.0 μl**  |                  |           |
 |**Add:**		       |             |                   ||
 | *Template DNA*	   |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |

|Thermal  |cycling  |conditions|
|:------|:---------|:----------|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***60°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   


**Primers: CytB_F1_Sequence (5' to 3'):** CTATTAATTTAGTTAAAGCACAC; **CytB_R1_Sequence (5' to 3'):** ACAGAATAATCTCTAGCACCA


**Table 6.4 *mdr1***

|*Pfmdr1 (4155bp)*	  | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                           
|:------------------|:--------------:|:-------------|:--------------:|
| *Master Mix:*      |                                                                                                 	                                               
| *5X GC Buffer*    | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)* |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H2O*	            |  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl** |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |

|Thermal  |cycling  |conditions|
|:----------:|:---------:|:----------:|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***60°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   


**Primers:  Pfmdr1_F1_Sequence (5' to 3'):** TGGTAACCTCAGTATCAAAG; **Pfmdr1_R1_Sequence** (5' to 3'): CATCTTGTGCTGATAATAATTC





**Table 6.5 *dhfr***

|*Pfdhfr (2067bp)*| **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                          
|:------------------:|:--------------:|:-------------:|:--------------:|
| *Master Mix:*      |                                                                                                 	                                      
| *5X GC Buffer*    | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)* |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H<sub>2</sub>O*	            |  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl**	 |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	        |**50.0 μl**	|	                 |          |

|Thermal  |cycling  |conditions|
|:------:|:---------:|:----------:|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***58°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   

**Primers:  Pfdhfr_F1 Sequence (5' to 3'):** TTTTTACTAGCCATTTTTGTATTCC; **Pfdhfr_R1 Sequence  (5' to 3'):** TTAACCGTTCAGGTAATTTTGTCA

***\*Primers adapted from: SC, Carlton JM. 2016.*** A Method for Amplicon Deep Sequencing of Drug Resistance Genes in Plasmodium falciparum Clinical Isolates from India. J Clin Microbiol 54:1500–1511.





**Table 6.6 *dhps***

|*Pfdhps (2817bp)*	  | **50.0 μl rxn**| X samples+1 |	***Final [conc]***  |                                          
|:------------------|:--------------|:-------------|:--------------|
| *Master Mix:*      |                                                                                                 	                                                  
| *5X GC Buffer*    | 	*10.0 μl*|                  | *1 x*|
| *dNTPs (10mM)*	  |   *1 μl*	 |                  | *0.2 mM*|
| *F primer (10uM)* |	*1.25 μl*	 |                  | *0.25 uM*|
| *R primer (10uM)*	|  *1.25 μl* |                  |*0.25 uM*|
 |*H<sub>2</sub>O*|  *32 μl*   |                  |         |
 |**Add last:** 	  |            |                  |         |
 |*HF Phusion Taq*  |	  *0.5 μl*  |                  |*1 unit*|
 |*Total*	          | **46.0 μl** 	 |                  |           |
 |**Add:**		      |            |                   ||
| *Template DNA*	  |*4.0 μl/well*| 	 ------	      |*50-250 ng*|
 |**TOTAL**	          |**50.0 μl**	 |	                 |          |
 
 |Thermal  |cycling  |conditions|
|:---------:|:---------:|:----------:|
|*98°C*	  |*3:00 min*	   |           |
|***98°C***	|***0:10 min***|           |
|***58°C***	|***0:30 min***	 | ***30X*** |
|***65°C***	|***5:00 min***	 |       |
|*65°C*	|*10:00min*	 |           |
|*4°C* |*∞*|        |   


**Primers:  Pfdhps_F Sequence (5' to 3'):** AATATTTGCGCCAAACTTTTTA; **Pfdhps_R Sequence (5' to 3'):** TTTATTTCGTAATAGTCCACTTTTGAT

***\*Primers adapted from: SC, Carlton JM. 2016.*** A Method for Amplicon Deep Sequencing of Drug Resistance Genes in Plasmodium falciparum Clinical Isolates from India. J Clin Microbiol 54:1500–1511.


<a id="chapter-6"></a>
## Electrophoresis 

This step is necessary to ensure successful amplification of amplicons. It is recommended to run at least 25% of the total samples, all no-template and negative controls on the gel to confirm amplification was successful and no contamination occured. Please note PCR amplification can be affected by numerous factors, including but not limited to, DNA quality and quantity. 

### Consumables
#### **Table7  Electrophoresis Consumables**

|Item                                                   |Quntity                |	Storage|
|:-----------------------------------------------------:|:-----------------------:|:------------------:|                                             
|**Agarose** |1g (for a 1% gel)|Room temperature|
|**1x solution of 10X TBE Buffer and deionized water** | 	100mL (for a 1% gel) | Room temperature|
**Nucleic Acid Gel Stain** |	5 μl per 100mL of buffer |	Room temperature
**Orange Dye** |	2 μl per 8 μl PCR product |	Room temperature

### Preparation
  - Create X solution of 10X TBE Buffer and deionized water before you start.

### Procedure   
1.  Choose an Erlenmeyer flask that is 2-4 times the volume of the solution and place a stirring rod into the flask.
2.  Weigh the agar to the desired concentration
      * For a 1% agarose gel, 1.0 gram of agarose + 100mL of buffer will fill a medium gel chamber
      * For a 1.5% agarose gel, 1.5 gram of agarose + 100mL of buffer will fill a medium gel chamber
3. 	Add the appropriate amount of buffer for the) settled agar.
6. 	Repeat steps 4-5 until all the agar is dissolved (no transparent agarose clumps should be present).
7. 	Allow the solution to cool on a stirring plate until you can comfortably hold the flask with your hands.
8. 	Using a 10 μl pipette, add nucleic acid gel stain to the solution. For every 100mL of buffer, add 5 μl of gel stain. Swirl solution to mix, making sure as little bubbles as possible are created.
9. 	Pour the cooled solution into the gel form- ensure no bubbles are present. Place the comb into the gel and allow the gel to sit undisturbed for at least 15 minutes or until the gel has become firm (the color will change from clear to slightly milky in color).
10. When gel has solidified, ensure the wells are aligned with the black (negative) nodes on the electrophoresis chamber and fill with buffer until it covers about a centimeter above the gel. Remove the comb.
11. Combine a mixture of 2 μl of orange dye and 8 μl of each sample and load 8uL of that mixture into each well.
12. Be sure to include reference ladders (no orange dye necessary).
13. Place the lid on the chamber box and connect the black node to the negative terminal and the red node to the positive terminal. Turn on the power supply and adjust the voltage to 100-130 volts.
14. Run gel for about 45 minutes -until the samples nearly reach the end of the gel.  **DO NOT** allow samples to run off the gel.
15. Turn off the power supply, disconnect the electrodes, and remove the lid.
16. Remove the gel from the chamber and take to the gel reading station for analysis.
17. Once amplification is confirmed, proceed to SequalPrep Normalization **(page 15-16)**.

<a id="chapter-7"></a>
## SequalPrep Normalization 
This step utilizes ChargeSwitch Technology to purify and normalize amplicon PCR product. PCR product is added to the plate and mixed with Binding Buffer, which then sits at room temperature for 1 hour. The DNA is then washed and eluted, resulting in a purified and normalized DNA product. 

This protocol format was adapted from the SequalPrep™ Normalization Plate (96) Kit [protocol](https://www.thermofisher.com/order/catalog/product/A1051001) from ThermoFisher Scientific

### Consumables


#### **Table 8.** PCR Purification 1 SequalPrep Consumables
|                Item                                                 | Quantity                |Storage                          |
|:------------------------------------------------------------------:|:------------------------:|:----------------------------------:|
|**SequalPrep™Normalization Plate (96)**    	                      |1 plate   	         | 15° to 30°C (Room temperature)  |
|**SequalPrep™Normalization Binding Buffer**                          |1 μl per 1 μl of sample   | 15° to 30°C (Room temperature) |
|**SequalPrep™Normalization Wash Buffer** 	                      |0 μl per sample 	         | 15° to 30°C (Room temperature) |
|**SequalPrep™Normalization Elution Buffer(10mM Tris-HCL, pH 8.5)**   |20 μl per sample          |15° to 30°C (Room temperature) |
|**(Optional) Foil Tape** 	                                      |2 sheets                  |15° to 30°C (Room temperature) |

### Preparation

 * You will need at least 250 ng amplicon per well to use with the SequalPrep™ Normalization Plate to achieve robust normalization. An average efficiency PCR (20ul reaction volume) produces in the range of 25-100ng/ul, allowing you to purify using 5-10ul using the SequalPrep system.  

 * Elution options: 
     - The **standard elution** method (described below) is designed to elute purified DNA from each well using 20 μl elution volume to obtain each amplicon at a concentration of 1–2 ng/μl. 
     - The **optional sequential** elution method is designed to sequentially elute multiple rows or columns using the same 20 μl of elution buffer to obtain higher amplicon concentrations. The amplicon concentrations will be additive as sequential wells are eluted. For example, dispense 20 μl of elution buffer into the first column (A1–H1), mix well, and incubate for 5 minutes at room temperature. Then, simply move this column of elution buffer to the next column (A2–H2), and again incubate for 5 minutes. Continue this step to obtain your specific elution needs for the downstream application of choice. 
     

   **NOTE:** Proceed very cautiously during this procedure and take your time- do not introduce bubbles when pipetting. 


### Procedure

#### **Binding Step**

1. 	Transfer the desired volume of PCR product (5–25 μl PCR reaction mix, at least 250 ng amplicon/well) from the PCR plate into the wells of the SequalPrep™ Normalization plate. 
2. 	Add an equivalent volume of SequalPrep™ Normalization Binding Buffer. 
    - For example: To purify 10 μl of PCR product, add 10 μl SequalPrep™ Normalization Binding Buffer.
3. 	Mix completely by pipetting up and down 10 times, or seal the plate with Foil Tape, vortex to mix, and briefly centrifuge the plate.
4. 	Incubate the plate for 1 hour at room temperature to allow binding of DNA to the plate surface. Mixing is not necessary at this stage.
       **NOTE:** Incubations longer than 60 minutes do not improve results. However, depending on your work flow you may perform overnight incubation at room temperature for the binding step.  

5. 	Optional: If >25 ng DNA/well yield is desired, transfer the amplicon/Binding Buffer mixture from Step 4 to another, fresh well/plate to sequentially bind more DNA. Perform DNA binding at room temperature for 1 hour. 
        **NOTE:** After binding is complete, you can remove the amplicon/Binding Buffer mixture from the well and store at –20ºC for up to 30 days to perform additional purifications at a later time.  
        
#### **Washing Step**

6. 	Aspirate the liquid from wells. Be sure not to scrape the well sides during aspiration.
7. 	Add 50 μl SequalPrep™ Normalization Wash Buffer to the wells. Mix by pipetting up and down twice to improve removal of contaminants.
        **NOTE:** If you wish to store the amplicon/Binding Buffer mixture for additional purifications at a later time, aspirate the liquid from wells into another plate and store at –20ºC for up to 30 days.  

8. 	Completely aspirate the buffer from wells and discard. To ensure complete removal of wash buffer and maximize elution efficiency, you may need to invert and tap the plate on paper towels depending on the pipetting technique or instrument used. A small amount of residual Wash Buffer (1–3 μl) is typical and does not affect the subsequent elution or downstream applications.

#### **Elution Step**

9. 	Add 20 μl SequalPrep™ Normalization Elution Buffer to each well of the plate.
       **NOTE:** Do not use water for elution. If you need to elute in any other buffer, be sure to use a buffer of pH 8.5–9.0. If the pH of the buffer is <8.5, the DNA will not elute efficiently.  

10. Mix by pipetting up and down 5 times (or seal the plate with Foil Tape, vortex to mix, and briefly centrifuge the plate). Ensure that the buffer contacts the entire plate coating (up to 20 μl level).
11. Incubate at room temperature for 5 minutes.
12. Transfer the normalized and purified DNA to a new 96-well skirted plate. You can store the eluted DNA at 4°C (short-term storage) or –20°C (long-term storage) until further use. 

13. Pool each drug resistance gene (normalized PCR amplicon) as follows: 
	  - ***10uL of each amplicon (i.e. mitochondria/cyt-b, k13, mdr1, dhfr, dhps, and crt) for each respective sample*** 

Expected DNA concentration is 1–2 ng/μl when using 20 μl elution volume. 



**SAFE STOPPING POINT** If you do not immediately proceed to *Tagment Genomic DNA,* seal plate with Microseal “B” adhesive seal and store it at ‐15° to ‐25°C for up to a week.

<a id="chapter-8"></a>
## Tagment Genomic DNA 
  This step uses the Nextera transposome to tagment gDNA, which fragments and then tags the DNA with adaptor sequences in a single step. 

![alt text](https://github.com/CDCgov/MaRS/blob/master/images/Nextera.png)




### Consumables

#### **Table 9.Tagment Genomic DNA Consumables**

|Item	|Quantity|	Storage|
|:-----------------------------------------:|:--------------------------:|:----------------:|
|**ATM (Amplicon Tagment Mix)**              |	5 μl per sample          |	-15° to -25°C|
|**TD (Tagment DNA Buffer)**	                |10 μl per sample           |	-15° to -25°C|
|**NT (Neutralize Tagment Buffer)**         	| 5 μl per sample       	    |15°C to 30°C|
|**TruSeq Index Plate Fixture (FC-130-1005)**|	1 (if available)            |	Room temperature|
|**96-well 0.2 ml PCR plate**                |	1 plate	                   |Room temperature|
|**Microseal “A” film**                     |	1	                          |Room temperature|


### Preparation
   - If not completed already, carefully pool 10uL of each gene’s corresponding sample into the same well of a half skirt plate before continuing
   - Be sure all samples are mixed thoroughly by pipetting gently 10 times.   
   - Thaw on ice, ATM and TD, invert thawed tubes 3-5 times and then centrifuge briefly. 
   - Check NT for precipitates. If present, vortex until all particulates are re - suspended. 
   - Set up thermal cycler and choose preheat lid option:   ***Thermocycler Program:*** ***55 °C for 5 min***, ***10 °C for ∞***


### Procedure
1. Add the following items in the order listed to each well of a new Hard-Shell skirted PCR plate. Pipette to mix.

|      Item                           |Volume (μl)|
|:------------------------------------:|:----------:|
|   **TD**	                          |   10  |
| **Normalized pooled gene amplicons**|   5   |

2.	Add 5 μl ATM to each well. Pipette to mix.
3.	Centrifuge at 280 × g at 20°C for 1 minute.
4.	Place on the preprogrammed thermal cycler and run the tagmentation program
	***Thermocycler Program:*** ***55 °C for 5 min***, ***10 °C for ∞***
	
  **NOTE:** Be sure to remove the plate as soon as the reaction has completed- NT must be added to the well immediately after the reaction has completed in order to stop the reaction.  

5.	Add 5 μl NT to each well. Pipette to mix.
6. 	Centrifuge at 280 × g at 20°C for 1 minute.
7.  Incubate at room temperature for 5 minutes.
8.	[optional] To assess tagmentation, run 1 uL on Agilent Bioanalyzer and/or Tapestation using high sensitivity DNA chip. 
9.	Proceed to Library Amplification and Index PCR

<a id="chapter-9"></a>
## Library Amplification and Index PCR 
This step amplifies the tagmented DNA using a limited-cycle PCR program. The PCR step adds Index 1 (i7) adapters and Index 2 (i5) adapters and sequences required for cluster formation. Use the full amount of recommended input DNA. To ensure libraries produce high-quality sequencing results, use the specified number of PCR cycles

### Consumables
#### **Table 10. Library Amplification and Index PCR Consumables**

|Item	                                                                                       |Quantity          |	Storage|
|:----------------------------------------------------------------------------------------------|:-----------------|:-----------------|
|**NPM (Nextera PCR Master Mix)**	                                                               |15 μl per sample  |	-15° to -25°C   |
|**Nextera XT Index 1 Primers (N7XX) from the Nextera XT Index kit (FC-131-1001 or FC-131-1002)**|5 μl per sample   |	-15° to -25°C   |
|**Nextera XT Index 2 Primers (S5XX) from the Nextera XT Index kit (FC-131-1001 or FC-131-1002)**|5 μl per sample	  |-15° to -25°C    |
|**96-well 0.2 ml PCR plate**                                                                   |1 plate           |	Room temperature|
|**Microseal “A” film**                                                                        	|1	              |Room temperature|


### Procedure
1. 	Using the entire 25μl from the Tagment Genomic DNA step (see page 17-18), prepare the following: .
2. 	Arrange the Index 1 and 2 primers in a rack (i.e. the TruSeq Index Plate Fixture) using the following arrangements as needed:
    - Arrange Index 2 primer tubes (white caps, clear solution) vertically, aligned with rows A through H.
    - Arrange Index 1 primer tubes (orange caps, yellow solution) horizontally, aligned with columns 1 through 12. For more information on index selection, see Dual Indexing Principle, on page 23.

*\* If no TruSeq Index Plate Fixture is available, arrange the Index adaptors in the same way, and individually pipette each adaptor into its corresponding well*

**Figure 2 *TrueSeq Index Plate Fixture***

![alt text](https://github.com/CDCgov/MaRS/blob/master/images/TruSeqIndex.png)             

         A Index 2 primers (white caps), B Index 1 primers (orange caps), C 96‐well plate.
                




3. 	Place the 96‐well PCR plate with the 25 μl of resuspended PCR product DNA in the TruSeq Index Plate Fixture. 
    -	Add 5 μl of each index into its corresponding well. 

4. 	Add 15uL NPM to each well containing index adaptors. Pipette to mix.  
5.	Centrifuge at 280 × g at 20°C for 1 minute.
6.	Cover the plate with Microseal 'A'.
7.  Perform PCR on a thermal cycler using the following program:
    * 72°C for 3 minutes
    * 95°C for 30 seconds
    * 12 cycles of:
        * 95°C for 10 seconds
        * 55°C for 30 seconds
        * 72°C for 30 seconds
    * 72°C for 5 minutes
    * Hold at 10°C
8.	If you choose to continue, proceed to PCR Clean-up 2




**SAFE STOPPING POINT** If you do not immediately proceed to Library Clean-Up on page 21, seal plate with an adhesive seal and store it at 2° to 8°C for up to a week.

<a id="chapter-10"></a>
## Library Clean‐Up 
This step uses AMPure XP beads to clean up the final library before quantification.

### Consumables

**Table 11. PCR Purification #2 Consumables**

|Item	                                  |Quantity	                 |Storage                                            |
|:--------------------------------------|:-------------------------|:--------------------------------------------------------|
|**RSB (Resuspension Buffer)**          |	52.5 μl per sample|-15° to -25°C (after initial thaw, can keep at 2° to 8°C)|
|**AMPure XP beads**                    |	90 μl per 50 μl of sample|	2° to 8°C|
|**Freshly Prepared 80% Ethanol (EtOH)**|	400 μl per sample|	Room temperature|
|**96‐well 0.2 ml PCR plate**           |	1 plate                  |                                                         |
|**[Optional] Microseal 'B' film**	    |	                         |                                                         |
|**96‐well MIDI plate**                 |1 plate|	                    |


### Preparation

• Determine whether or not a plate transfer is necessary.  If the PCR reaction volume multiplied by 2.8 exceeds the volume of the PCR plate, a transfer to a 300 μL round bottom plate or a 1.2 mL deep-well plate is required. 

• Thaw RSB at room temperature.  

• Bring the AMPure XP beads to room temperature- wait at least 30 minutes.

• Once at room temperature, shake the Agencourt AMPure XP bottle to re-suspend any magnetic particles that may have settled. Ensure magnetic beads are well (evenly) distributed before adding them to samples. 

**NOTE:** Proceed very cautiously during this procedure and take your time. 

**NOTE:** 70% ethanol is hygroscopic. That is, when opened the ethanol will both evaporate and absorb water over time. Re-use eventually will be at a lower concentration. There is also miscibility involved with ethanol and water. For example, measuring out 70 mL of ethanol and topping off to 100 mL with water will generate ~65% ethanol. Measuring 80 mL ethanol and 20 mL water separately, then combining them will generate ~95 mL of 80% ethanol. Make sure to use molecular biology grade water (DNAse, RNase and Protease free). 

### Procedure
1. 	Centrifuge the Library Amplification plate at 1,000 × g at 20°C for 1 minute to collect condensation, carefully remove seal.
2. 	Using a multichannel pipette, transfer the entire 50 μl of PCR product from the Library Amplification plate to the MIDI plate. Change tips between samples.

**NOTE:** Transfer the sample to a 96‐well MIDI plate if planning to use a shaker for mixing. If mixing by pipette, the sample can remain in the 96‐well PCR plate.  

3. 	Gently shake the AMPure XP beads for 30 seconds to make sure that the beads are evenly dispersed. Add an appropriate volume of beads to a trough depending on the number of samples being processed and desired fragment selection. Smaller amplicons in Nextera XT library preps typically yield smaller insert size ranges. 
To maximize recovery of smaller fragments from the bead cleanup step, use the following conditions: 

![alt text](https://github.com/CDCgov/MaRS/blob/master/images/AMPure.png)

4.  Using a multichannel pipette, add 90 μl of AMPure XP beads to each well of the Amplicon PCR plate. Change tips between columns. 
5.  Gently pipette entire volume up and down 10 times if using a 96‐well PCR plate, or seal plate and shake at 1800 for 2 minutes if using a MIDI plate. Change tips after each column. The mixture should appear homogeneous. 
6.  Incubate the mixed samples at room temperature without shaking for 5 minutes.
7.  Place the Library amplification plate on a magnetic stand for 2 minutes. WAIT for the solution to clear before proceeding.
8.  With the Amplicon PCR plate on the magnetic stand, use a multichannel pipette to carefully remove and discard all the supernatant. Change tips between samples.

	**DO NOT disturb the ring of separated magnetic beads.** 
	
9.  With the Library amplification plate on the magnetic stand, wash the beads with freshly prepared 80% ethanol as follows:
      * Using a multichannel pipette, add 200 μl of freshly prepared 80% ethanol to each sample well.
      * Incubate the plate on the magnetic stand for 30 seconds at room temperature.
      * Carefully remove and discard the ethanol. 
    
     **Note: The beads are not drawn out easily when in alcohol, so it is not necessary to leave any supernatant behind.**

10. With the Library amplification plate on the magnetic stand, perform a second ethanol wash as follows:
      * Using a multichannel pipette, add 180 μl of freshly prepared 80% ethanol to each sample well.
      * Incubate the plate on the magnetic stand for 30 seconds.
      * Carefully remove and discard all the ethanol.
      * Use a P20 multichannel pipette with fine pipette tips to remove excess ethanol.
11. With the Library amplification plate still on the magnetic stand, allow the beads to air‐dry for 15 minutes.
    
    **NOTE: make sure not to over dry the beads. Bead pellets will appear cracked if over dried.** 
  
12.    Remove the Amplicon PCR plate from the magnetic stand. Using a multichannel pipette, add 52.5 μl RSB to each well of the Amplicon PCR plate.
13.    Gently pipette mix up and down 10 times, or seal plate and shake at 1800 for 2 minutes if using a MIDI plate. Change tips after each column. 
14.    Incubate at room temperature for 2 minutes.
15.    Place the plate back on the magnetic stand for 2 minutes or until the supernatant has cleared.
16.    Using a multichannel pipette, carefully transfer 50 μl of the supernatant from the Library amplification plate to a new 96‐well PCR plate. Change tips between samples to avoid cross‐contamination.


**SAFE STOPPING POINT**
If you do not plan to proceed to Library Clustering on page 24, seal the plate with Microseal “B” adhesive seal. Store the plate at ‐15° to ‐25°C for up to a week.

<a id="chapter-11"></a>
## Library Clustering 
It is important to consider library size when preparing samples for cluster generation. Because the clustering process preferentially amplifies shorter libraries in a mixture of fragments, large libraries tend to cluster less efficiently than small libraries. The DNA concentration used for clustering can be adjusted to increase the cluster density of larger libraries.  Consider table 1 below:

### Library Denaturing and MiSeq Sample Loading

##### Guidelines for Optimal Cluster Density

![Cluster](https://github.com/CDCgov/MaRS/blob/master/images/ClusterDensity.png)

**Average Library Size** |  **Conversion Factor** |	**DNA Concentration for Cluster Generation**
-----------| ----------------| -----------------
250 bp |	1 ng/μl = 6 nM |	6-12 pM
500 bp | 	1 ng/μl = 3 nM | 	6-12 pM
1,000 – 1,500 bp |	1 ng/μl = 1.5 nM |	12-20 pM

*The values presented here are approximations, and exact values determined for each experiment may differ from these guidelines. The guidelines presented are applicable to Nextera DNA libraries and Nextera XT libraries that have not been normalized. 

<a id="chapter-12"></a>
## Library Pooling, Quantification, and Normalization 
#### This step requires three parts: 
   **Part I** Pool libraries
   **Part II** Quantification of fragment size and concentration to determine library concentration in nM
   **Part III** Diluting your final library in Resuspension Buffer (RSB) or fresh 10 mM Tris pH 8.5 to a 4 nM solution
     
### **Part I:** Pooling
Aliquot 5 μl of diluted DNA from each library into a 1.5 microcentrifuge tube and mix aliquots for pooling libraries with unique indices. Depending on coverage needs, up to 384 libraries can be pooled for one MiSeq run.

### **Part II:** Quantification
Illumina recommends quantifying your libraries using a fluorometric quantification method that uses dsDNA binding dyes.

* In order to determine the fragment size, this laboratory adopted the Agilent D5000 ScreenTape System Quick Guide protocol from Agilent Technologies. 

* In order to determine the library concentration, this laboratory adopted the Qubit® dsDNA HS Assay Kits protocol from Life Technologies. 

### DNA Concentration in nM
After determining the fragment size and concentration of your pooled product, you will calculate the DNA concentration in nM, based on the size of DNA amplicons as determined by an Agilent Technologies 2100 Bioanalyzer trace and concentration by Qubit:

```
(concentration in ng/μl) * (10^6) / (660g/mol) * (average library size) = concentration in nM
```

_For example:_

(15ng/μl) x (10^6) / (660g/mol) x (500bp) = 45 nM

### **Agilent Technologies Agilent D5000 ScreenTape System** 
This SOP format was adapted from the Agilent D5000 ScreenTape System Quick Guide protocol from Agilent Technologies. 

#### Consumables
###### **Table 12. TapeStation Consumables**
Item |  Quantity |	Storage
-----| --------| -------
**Sample Buffer** |	10 μl per sample |	2°– 8° C
**D5000 Ladder** | 	1 μl | 2°– 8° C
**ScreenTape** |	Holds 16 samples per tape |	2°– 8° C

#### Procedure

##### Prepare TapeStation System D5000

* Launch the 2200 TapeStation Controller Software.
* Load single D5000 ScreenTape device and loading tips into the 2200 TapeStation instrument.

##### Sample Preparation D5000 ScreenTape Assay   
1.	Allow reagents to equilibrate at room temperature for 30 minutes. 
2.	Vortex mix before use.
3.	Prepare ladder by mixing 10 μl D5000 Sample Buffer (green lid) with 1 μl D5000 Ladder (yellow lid) in a tube strip.
4.	Prepare sample by mixing 10 μl D5000 Sample Buffer (green lid) with 1 μl DNA sample in different tube strips.
5.	Spin down, then vortex using IKA vortexer and adaptor at 2000 rpm for 1 minute.
6.	Spin down to position the sample at the bottom of the tube.

##### Sample Analysis
1.	Load samples into the 2200 TapeStation instrument
2.	Select the required samples on the 2200 TapeStation Controller Software.
3.	Click Start and specify a filename with which to save your results. 


##### SAFE STOPPING POINT
If you do not plan to proceed to Part II Qubit Flurometer 3.0 dsDNA HS Assay on page 26, leave your sample in 4°C. 

### Qubit Fluorometer 3.0 dsDNA HS Assay 
This SOP format was adapted from the Qubit® dsDNA HS Assay Kits protocol from Life Technologies. 

#### Consumables
###### **Table 13. Qubit 3.0 Fluorometer Consumables**
Item |  Quantity |	Storage
--------------| ---------------------| -------------
**Qubit dsDNA HS Buffer** |	199 μl per sample for working solution |	Room temperature
**Qubit dsDNA HS Reagent** | 	1 μl per 199 μl of HS Buffer | Room temperature
**Standard #1** |	10 μl per use  |	 2°– 8° C
**Standard #2** |	10 μl per use  |	 2°– 8° C
**Qubit™ Assay Tubes** |	1 per sample and 1 for each |	Room temperature

#### Before you begin 
* The final volume in each tube must be 200 μl. 
* Each standard tube requires 190 μl of Qubit working solution + 10 μl of the standard 
* Each sample tube requires anywhere from 180–199 μl + the corresponding volume to complete the necessary 200 μl: this laboratory uses 195 μl working solution + 5 μl of sample
* Careful pipetting is critical to ensure that the exact volume is added to the working solution—work SLOWLY
* Be sure to use a clean plastic tube each time you prepare Qubit working solution.  Do not mix the working solution in a glass container 

#### Procedure
##### Standard and Sample Preparation
1. Prepare the tubes:

     * Set up two (2) 0.5-mL tubes for standards, and the required number of tubes for samples. 
          * **Note** Use only the thin-wall, clear, 0.5-mL PCR tubes (described in Table 2 User‐Supplied Consumables)

     * Label the tube lids- do not label the side of the tube as this could interfere with the sample read

2. Prepare the Qubit working solution:

     * Prepare sufficient Qubit working solution to accommodate all standards and samples by diluting the Qubit dsDNA HS Reagent 1:200 in Qubit dsDNA HS Buffer.

          * 1 μl Qubit dsDNA HS Reagent + 199 μl Qubit dsDNA HS Buffer
     
          * For example, for 8 samples, prepare enough working solution for the samples and 2 standards: ~200 μl per tube in 10 tubes yields 2 mL of working solution (10 μl of Qubit reagent plus 190 μl of Qubit buffer).
     
3. Prepare the standards:

     * Add 190 μl of Qubit working solution to each of the tubes used for standards. 

     * Add 10 μl of each Qubit standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.

4. Prepare the samples:

     * Add Qubit working solution to individual assay tubes so that the final volume in each tube after adding the sample is 200 μl.
     * Add each sample to the assay tubes containing the correct volume of Qubit working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μl.
     
        * **Note** Your sample can be anywhere from 1–20 μl. Add a corresponding volume of Qubit working solution to each assay tube: anywhere from 180–199 μl.

| |   |   |   |   |  |
|------ | ---- | ---- | ---- | ---- | ---- |
|Working Sample Volume | 199 μl | 195 μl | 190 μl | 185 μl | 180 μl|
|Sample Volume  | 1 μl | 5 μl | 10 μl | 15 μl | 20 μl |

5. Allow all tubes to incubate at room temperature for 2 minutes. 

##### Standard and Sample Reading

1. On the home screen of the Qubit 3.0 Fluorometer, select “**dsDNA**”, then “**High Sensitivity**”, and then “**Read Standards**.”  
     * **Note**: If you have already performed a calibration for the selected assay, the instrument prompts you to choose between reading new standards and running samples using the previous calibration. If you wish to use the previous calibration, disregard step 3 in the Standard and Sample Preparation step, and skip to step 4 below. Otherwise, continue with step 2 below. 
          
2. Insert the tube containing Standard #1 into the sample chamber, close the lid, and then press **Read standard**. When the reading is complete (~3 seconds), remove Standard #1. 
     
3. Insert the tube containing Standard #2 into the sample chamber, close the lid, and then press **Read standard**. When the reading is complete, remove Standard #2.  
     
4. Press **Run samples**.
     
5. On the assay screen, select the sample volume and units using the + or – buttons on the wheel to select the sample volume added to the assay tube (from 1–20 μl).
     
6. From the dropdown menu, select the units for the output sample concentration (ng/ μl). 
     
7. Insert a sample tube into the sample chamber, close the lid, and press Read tube. When the reading is complete (~3 seconds), remove the sample tube and repeat until all samples have been read. 

### **Part III:** Normalization
Dilute concentrated final library using Resuspension Buffer (RSB) or fresh 10 mM Tris pH 8.5 to **4 nM**. 

Example:
Given a calculated concentration of 45nM, use C<sub>1</sub>V<sub>1</sub>=C<sub>2</sub>V<sub>2</sub> to calculate how much RSB and sample to mix to create a 4nM concentration:

45nM (V<sub>1</sub>)= 4nM (20 μl)

V<sub>1</sub>  = 1.78 μl of sample + 18.22 μl of RSB produces 20 μl of a 4nM concentration

##### SAFE STOPPING POINT
If you do not plan to proceed to Library Denaturing and MiSeq Sample Loading on page 27, leave your sample in 4°C. 

<a id="chapter-13"></a>
## Library Denaturing and MiSeq Sample Loading 
In preparation for cluster generation and sequencing, pooled libraries are denatured with NaOH, diluted with hybridization buffer, and then heat denatured before MiSeq sequencing. Each run must include a minimum of 5% PhiX to serve as an internal control for these low-diversity libraries. Illumina recommends using MiSeq v2 reagent kits for improved run metrics.

#### Consumables
###### **Table 14.** Library Denaturing and MiSeq Sample Loading Consumables
Item |  Quantity |	Storage
-----| --------| -------
**RSB (Resuspension Buffer)** |	6 μl |	-15° to -25°C
**HT1 (Hybridization Buffer)** | 	1540 μl | -15° to -25°C
**0.2 N NaOH (less than a week old)** |	10 μl |	Room temperature
**200mM Tris-HCl pH7.0** | 5 μl | Room temperature
**PhiX Control Kit v3 (FC‐110‐3001)** | 2 μl | -15° to -25°C
**MiSeq v2 reagent cartridge** | 1 cartridge | -15° to -25°C
**1.7 ml microcentrifuge tubes** | 3 tubes | Room temperature
**2.5 L ice bucket** | 1 bucket | Room temperature

#### Procedure

##### Preparation 
1. Begin thawing the reagent cartridge and HT1 before denaturing and diluting libraries by placing them in a room temperature water bath for about an hour 
     * Once thawed, store the cartridge and HT1 in the ice bucket until ready for sample loading. 
2. Obtain an ice bucket for your thawed cartridge, freshly made reagents, and sample.
3. Check pH of the stock 1.0N NaOH and the resulting 0.2N NaOH dilution using pH reader. 
     * **Note** CO2 in the room will acidify NaOH over time. It is absolutely **critical** that the NaOH has a minimum pH >12.5. 
4. Prepare a fresh dilution of 0.2 N NaOH [this is a critical step; NaOH must be prepared fresh every time] 
     * Using a 1000ul pipette, measure out 800uL of laboratory-grade water.
     * In a separate microcentrifuge tube, measure 200uL of stock 1.0N NaOH.
     * Combine the two volumes and then invert several times to mix
          * **Note** This results in a 1mL of 0.2N NaOH solution; the resulting solution cannot be stored and must be used within 6 hours
          * **Note** The final concentration of NaOH cannot exceed 0.001 (1mM). Higher NaOH concentrations will inhibit library hybridization to the flow cell and result in very low cluster density. 
5. If you have not already done so, prepare a 200mM stock of Tris-HCl pH7.0 by combining 800 μl of Laboratory-grade water and 200 μl of Tris-HCl 1M.

##### Denature DNA

1. Combine the following volumes of pooled final DNA library and freshly diluted 0.2 N NaOH in a microcentrifuge tube:
     * 4 nM pooled library (5 μl)
     * 0.2 N NaOH (5 μl) 
2. Set aside the remaining dilution of 0.2 N NaOH to prepare a PhiX control within the next 12 hours.
3. Vortex briefly to mix the sample solution, and then centrifuge the sample solution at 280 × g (or about 1500rpm) at 20°C for 1 minute.
4. Incubate for 5 minutes at room temperature to denature the DNA into single strands.
5. To the 10 μl of denatured library, add 5 μl of 200mM Tris-HCl pH7.0 to neutralize the NaOH. 
6. Add the following volume of pre‐chilled HT1 to the tube containing denatured DNA:
     * Denatured DNA + Tris-HCl (15 μl)
     * Pre‐chilled HT1 (985 μl)
          * Adding the HT1 results in a 20 pM denatured library in 1 mM NaOH.
7. Place the denatured DNA on ice until you are ready to proceed to final dilution.

### Quick Review/Guide for denaturing 4nM library:

           5 μl of 4nM library + 5 μl of 0.2N NaOH + 5 μl of 200mM Tris-HCl pH 7.0

                                       Equals

                 15 μl of 1.3nM library and 0.067N NaOH + 66.7nM Tris pH 7.0

                              (Add 985 μl of chilled HT1)

                                       Equals

           1mL of 0.001N NaOH and 20pM denatured library + 1mM Tris-HCl pH 7.0
           
**NOTE: If you have to start with a lower concentration library, follow the below protocol for denaturing a 2nM library.**

### Quick Review/Guide for denaturing 2nM library:
             5 μl of 2nM library + 5 μl of 0.2N NaOH + 5 μl of 200mM Tris-HCl pH7.0 

                                       Equals

                15 μl of 0.67nM library and 0.067N NaOH + 66.7nM Tris pH 7.0

                             (Add 985 μl of chilled HT1)

                                       Equals

              1mL of 0.0005N NaOH and 10pM denatured library + 1mM Tris-HCl pH 7.0

#### **Dilution chart for 10pM library:**
Desired Final Concentration | 6pM | 8pM | 10pM 
------ | ---- | ---- | ---- 
**10pM denatured library** | 360 μl | 480 μl | 600 μl 
**Pre-chilled HT1** | 240 μl | 120 μl | 0 μl 

#### Dilute Denatured DNA

1. Dilute the denatured DNA to the desired concentration using the following table. Illumina recommends targeting 1000–1200 K/mm² raw cluster densities using MiSeq v2 reagents. 

#### **Final Library Dilution chart**
Final Concentration | 2pM | 4pM | 6pM | 8pM | 10pM | 12pM | 15pM | 20pM
------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
**20pM denatured library** | 60 μl | 120 μl | 180 μl | 240 μl | 300 μl | 360 μl | 450 μl | 600 μl
**Pre-chilled HT1** | 540 μl | 480 μl | 420 μl | 360 μl | 300 μl | 240 μl | 150 μl | 0 μl 

2. Invert several times to mix and then pulse centrifuge the DNA solution.
3. Place the denatured and diluted DNA on ice.

### Denature and Dilution of PhiX Control
Use the following instructions to denature and dilute the 10 nM PhiX library to the same loading concentration as the Amplicon library. The final library mixture must contain at least 5% PhiX.

1. Combine the following volumes to dilute the PhiX library to 4 nM:
     * 10 nM PhiX library (2 μl)
     * RSB (3 μl)
2. Combine the following volumes of 4 nM PhiX and 0.2 N NaOH in a microcentrifuge tube:
     * 4 nM PhiX library (5 μl)
     * 0.2 N NaOH (5 μl)
3. Vortex briefly to mix the 2 nM PhiX library solution.
4. Incubate for 5 minutes at room temperature to denature the PhiX library into single strands.
5. To the 10 μl of denatured library, add 5 μl of 200mM Tris-HCl pH7.0 to neutralize the NaOH. 
6. Add the following volumes of pre‐chilled HT1 to the tube containing denatured PhiX library to result in a 20 pM PhiX library:
     * Denatured PhiX library (15 μl)
     * Pre‐chilled HT1 (985 μl)
7. Dilute the denatured 20 pM PhiX library to the same loading concentration as the Amplicon library as follows: 

#### **PhiX Dilution chart**
Final Concentration | 2pM | 4pM | 6pM | 8pM | 10pM | 12pM | 15pM | 20pM
------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ----
**20pM denatured library** | 60 μl | 120 μl | 180 μl | 240 μl | 300 μl | 360 μl | 450 μl | 600 μl
**Pre-chilled HT1** | 540 μl | 480 μl | 420 μl | 360 μl | 300 μl | 240 μl | 150 μl | 0 μl 

8. Invert several times to mix and then pulse centrifuge the DNA solution.
9. Place the denatured and diluted PhiX on ice.

### Combine Amplicon Library and PhiX Control
The recommended PhiX control spike‐in of ≥ 5% for low diversity libraries is possible with RTA v1.17.28 or later, which is bundled with MCS v2.2. For optimal performance, update to v3 software (MCS 2.3). If you are using an older version of the MiSeq software or sequencing these libraries on the GA or HiSeq, Illumina recommends using ≥ 25% PhiX control spike‐in.

1. Combine the following volumes of denatured PhiX control library and your denatured amplicon library in a microcentrifuge tube, which results in a 5% spike-in of PhiX:
     * Denatured and diluted PhiX control (30 μl)
     * Denatured and diluted amplicon library (570 μl)
2. Set the combined sample library and PhiX control aside on ice until you are ready to load the mixture into the MiSeq v2 reagent cartridge.
3. Invert the tube 1–2 times to mix and load all 600ul into the designated well in the cartridge.

<a id="chapter-14"></a>
## Supporting Information 

The protocols described in this guide assume that you are familiar with the contents of this section and have obtained all of the requisite equipment and consumables.

### Acronyms
#### Table 15. Definitions and Acronyms
Acronym | Definition 
------ | ------------------------------ 
**PCR** | Polymerase Chain Reaction- a technique used to amplify 1 to a few copies of a piece of DNA across several orders of magnitude, generating thousands to millions of copies of a single DNA strand  
**Primer** | A strand of short nucleic acid sequences that serves as a starting point for DNA synthesis during PCR 
**Amplicon** | A piece of amplified DNA that is the product of a PCR reaction

### Dual Indexing Principle
The dual indexing strategy uses two 8 base indices, Index 1 (i7) adjacent to the P7 sequence, and Index 2 (i5) adjacent to the P5 sequence. Dual indexing is enabled by adding a unique Index 1 (i7) and Index 2 (i5) to each sample. The 96 sample Nextera XT Index Kit (FC‐131–1002) use 12 different Index 1 (i7) adapters (N701–N712) and 8 different Index 2 (i5) adapters (S501–S508). The 24 sample Nextera XT Index Kit (FC‐131–1001) uses 6 different Index 1 (i7) adapters (N701–N706) and 4 different Index 2 (i5) adapters (S501–S504). In the Index adapter name, the N or S refers to Nextera XT sample preparation, 7 or 5 refers to Index 1 (i7) or Index 2 (i5), respectively. The 01–12 refers to the Index number. A list of index sequences is provided for generating sample sheets to demultiplex the samples

#### Index Pairs
Index 1 (i7) | Sequence | Index 2 (i5) | Sequence 
------ | ------ | ------ | ------ 
N701 | TAAGGCGA | S501 | TAGATCGC
N702 | CGTACTAG | S502 | CTCTCTAT
N703 | AGGCAGAA | S503 | TATCCTCT
N704 | TCCTGAGC | S504 | AGAGTAGA
N705 | GGACTCCT | S505 | GTAAGGAG
N706 | TAGGCATG | S506 | ACTGCATA
N707 | CTCTCTAC | S507 | AAGGAGTA
N708 | CAGAGAGG | S508 | CTAAGCCT
N709 | GCTACGCT |      |         
N710 | CGAGGCTG |      |         
N711 | AAGAGGCA |      |         
N712 | GTAGAGGA |      |         


### Low Plexity Pooling Guidelines
Illumina uses a green laser or LED to sequence G/T and a red laser or LED to sequence A/C. At each cycle, at least one of two nucleotides for each color channel are read to ensure proper registration. It is important to maintain color balance for each base of the index read being sequenced, otherwise index read sequencing could fail due to registration failure. If you choose the dual‐indexed sequencing workflow, always use at least two unique and compatible barcodes for each index (index 1 and index 2). The following tables illustrate possible pooling strategies:

#### **Table 16. Libraries Pooled: 6 or fewer; Sequencing Workflow; Single Index**
Plex |  Index 1 (i7) Selection |	Index 2 (i5) Selection
----------| ------------------------------| -------
**1-plex (no pooling)** |	Any index 1 adapter |	Any index 2 adapter
**2-plex** | 	[option 1]: N701 and N702;  [option 2]: N702 and N704
**3-plex** |	[option 1]: N701, N702 and N704; [option 2]: N703, N705, and N706
**4- or 5-plex** | [option 1]: N701, N702, N704, and any other Index 1 adapter; [option 2]: N703, N705, N706, and any other Index 1 adapter 
**6-plex** | N701, N702, N703, N704, N705, and N706 


#### **Table 17. Sequencing Workflow: Single or Dual Index**
Plex |  Index 1 (i7) Selection |	Index 2 (i5) Selection
---------------| -----------------------| ------------------
**7-12 plex, Dual Index** |	[option 1] N701, N702, N704, and any other Index 1 adapter (as needed); [option 2] N703, N705, N706, and any other Index 1 adapter (as needed) | [option 1] S501 and S502; [option 2] S503 and S504; [option 3] S505 and S506
**7-12 plex Single Index (96 sample Nextera Index Adapter kit)** | N701 - N706 and any other Index 1 adapter (as needed) | Any Index 2 (i5) adapter
**Greater than 12-plex** |	N701, N702, N703, N704, N705, N706, and any other Index 1 adapter | [option 1] S501, S502, and any other Index 2 adapter (as needed); [option 2] S503, S504, and any other Index 2 adapter (as needed); [option 3] S505, S506, and any other Index 2 adapter (as needed)

These strategies represent only some of the acceptable combinations. Alternatively, check the real sequences of each index in the tables to make sure that each base position has a signal in both color channels for the index read:


| Good Index 1 | Good Index 2 | Bad Index 1 | Bad Index 2 |
| ------: | ------: | ------: | ------: |
| 705  GG**A**CT**CC**T | 503  T**A**T**CC**T**C**T | 705  GG**AC**T**CC**T | 502  **C**T**C**T**C**T**A**T |
| 706  T**A**GG**CA**TG | 503  T**A**T**CC**T**C**T | 706  T**A**GG**CA**TG | 502  **C**T**C**T**C**T**A**T |
| 701  T**AA**GG**C**G**A** | 504  **A**G**A**GT**A**G**A** | 701  T**AA**GG**C**G**A** | 503  T**A**T**CC**T**C**T |
| 702  **C**GT**AC**T**A**G | 504  **A**G**A**GT**A**G**A** | 702  **C**GT**AC**T**A**G | 503  T**A**T**CC**T**C**T |
|      √√√√√√√√ |       √√√√√√√√ |      √√√√√√√√ |        √√√√xxxx |

√=signal in both color  
x=signal missing in one color channel  
(**bolded base** = color red)  
(non-bolded base = color green)  

### Prevent PCR Product Contamination
The PCR process is commonly used in the laboratory to amplify specific DNA sequences. Unless proper laboratory hygiene is used, PCR products can contaminate reagents, instrumentation, and genomic DNA samples, causing inaccurate and unreliable results. PCR product contamination can shut down lab processes and significantly delay normal operations.

Make sure that the lab is set up appropriately to reduce the risk of PCR product contamination:

* Physically Separate Pre-PCR and Post-PCR Areas
     * Physically separate laboratory space where pre‐PCR processes are performed (DNA extraction, quantification, and normalization) from the laboratory space where PCR products are made and processed (post‐PCR processes).
     * Never use the same sink to wash pre‐PCR and post‐PCR troughs.
     * Never share water purification systems for pre‐PCR and post‐PCR processes.
     * Store all supplies used in the protocols in the pre‐PCR area, and transfer to the post‐ PCR area as needed.
* Use Dedicated Equipment and Supplies
     * Dedicate separate full sets of equipment and supplies (pipettes, centrifuges, oven, heat block, etc.) to pre‐PCR and post‐PCR lab processes, and never share between processes.
     * Dedicate separate storage areas (freezers and refrigerators) to pre‐PCR and post‐PCR consumables.

Because the pre‐ and post‐amplification reagents are shipped together, it is important to unpack the reagents in the pre‐PCR lab area. After unpacking the reagents, move the post-amplification reagents to the proper post‐PCR storage area.

### Pre‐PCR and Post‐PCR Lab Procedures
To prevent PCR product contamination, it is important to establish lab procedures and follow best practices. Illumina recommends daily and weekly cleaning of lab areas using 0.5% Sodium Hypochlorite (10% Bleach).

###### CAUTION
To prevent sample or reagent degradation, make sure that all vapors from the cleaning solution have fully dissipated before beginning any processes.

### Daily Cleaning of Pre‐PCR Area
A daily cleaning of the pre‐PCR area using a 0.5% Sodium Hypochlorite (10% Bleach) solution helps to eliminate PCR product that has entered the pre‐PCR area. Identify pre‐PCR areas that pose the highest risk of contamination, and clean these areas with a 0.5% Sodium Hypochlorite (10% Bleach) solution before beginning any pre‐PCR processes. 
High‐risk areas might include, but are not limited to, the following items:

* Benchtops
* Door handles
* Refrigerator/freezer door handles
* Computer mouse
* Keyboards

### Daily Cleaning of Post‐PCR Area
Reducing the amount of PCR product in the post‐PCR area helps reduce the risk of contamination in the pre‐PCR area. Daily cleaning of the post‐PCR area using a 0.5% Sodium Hypochlorite (10% Bleach) solution helps reduce the risk of contamination. Identify post‐PCR areas that pose the highest risk of contamination, and clean these areas with a 0.5% Sodium Hypochlorite (10% Bleach) solution daily. 
High‐risk areas might include, but are not limited to, the following items:

* Thermal cyclers
* Bench space used to process amplified DNA
* Door handles
* Refrigerator/freezer door handles
* Computer mouse
* Keyboards

### Weekly Cleaning of All Lab Areas
One time a week, perform a thorough cleaning of the pre‐PCR and post‐PCR areas using 0.5% Sodium Hypochlorite (10% Bleach).

* Clean all benchtops and laboratory surfaces.
* Clean all instruments that are not cleaned daily.
* Thoroughly mop lab floors.
* Make sure that personnel responsible for weekly cleaning are properly trained on prevention of PCR product contamination.

###  Items Fallen to the Floor
The floor is contaminated with PCR product transferred on the shoes of individuals coming from the post‐PCR area; therefore, anything falling to the floor must be treated as contaminated.

* Disposable items that have fallen to the floor, such as empty tubes, pipette tips, gloves, lab coat hangers, must be discarded.
* Non‐disposable items that have fallen to the floor, such as a pipette or an important sample container, must be immediately and thoroughly cleaned. Use a 0.5% Sodium Hypochlorite (10% Bleach) solution to remove PCR product contamination.
* Clean any lab surface that has come in contact with the contaminated item. Individuals handling anything that has fallen to the floor, disposable or non‐disposable, must discard their lab gloves and put on a new pair.

### Best Practices
When preparing libraries for sequencing, always adhere to good molecular biology practices. Read through the entire protocol before starting to make sure that all of the required materials are available and your equipment is programmed and ready to use.

### Handling Liquids
Good liquid handling measures are essential, particularly when quantifying libraries or diluting concentrated libraries for making clusters.

* Small differences in volumes (±0.5 μl) can sometimes cause large differences in cluster numbers (~100,000).
* Small volume pipetting can be a source of potential error in protocols requiring the generation of standard curves, such as qPCR, or small but precise volumes, such as the Agilent Bioanalyzer.
* If small volumes are unavoidable, use due diligence to make sure that pipettes are correctly calibrated.
* Make sure that pipettes are not used at the volume extremes of their performance specifications. 
* Prepare the reagents for multiple samples simultaneously, to minimize pipetting errors, especially with small volume enzyme additions. As a result, pipette one time from the reagent tubes with a larger volume, rather than many times with small volumes. Aliquot to individual samples in a single pipetting movement to allow for standardization across multiple samples.

### Handling Magnetic Beads
NOTE Cleanup procedures have only been validated using the 96‐well plates and the magnetic stand specified in Tables 1 and 2.

Comparable performance is not guaranteed when using a microcentrifuge tube or other formats, or other magnets.

* Before use, allow the beads to come to room temperature.
* Do not reuse beads. Always add fresh beads when performing these procedures.
* Immediately before use, vortex the beads until they are well dispersed and the color of the liquid is homogeneous.
* When pipetting beads, pipette slowly and dispense slowly due to the viscosity of the solution.
* Take care to minimize bead loss, which can affect final yields.
* Change the tips for each sample, unless specified otherwise.
* Let the mixed samples incubate at room temperature for the time indicated in the protocol for maximum recovery.
* When removing and discarding supernatant from the wells, use a single channel or multichannel pipette and take care not to disturb the beads
* When aspirating the cleared solution from the reaction plate and wash step, it is important to keep the plate on the magnetic stand and not disturb the separated magnetic beads. Aspirate slowly to prevent the beads from sliding down the sides of the wells and into the pipette tips.
* To prevent the carryover of beads after elution, approximately 2.5 μl of supernatant is left when the eluates are removed from the bead pellet.
* Be sure to remove all of the ethanol from the bottom of the wells, as it can contain residual contaminants.
* Keep the reaction plate on the magnetic stand and let it air‐dry at room temperature to prevent potential bead loss due to electrostatic forces. Allow for the complete evaporation of residual ethanol, because the presence of ethanol affects the performance of the subsequent reactions. Illumina recommends at least minutes drying time, but a longer drying time can be required. Remaining ethanol can be removed with a 10 μl pipette.
* Avoid over drying the beads, which can impact final yields.
* Do not scrape the beads from the edge of the well using the pipette tip.
* To maximize sample recovery during elution, incubate the sample/bead mix for 2 minutes at room temperature before placing the samples onto the magnet.

### Avoiding Cross‐Contamination
Practice the following to avoid cross‐contamination:

* Open only one adapter tube at a time.
* Change the tips for each sample, unless specified otherwise.
* Pipette carefully to avoid spillage.
* Clean pipettes and change gloves between handling different adapter stocks.
* Clean work surfaces thoroughly before and after the procedure.

### Potential DNA Contaminants

When handling and processing samples using this protocol, use best practices to avoid PCR contamination, as you would when preparing PCR amplicons.

### Temperature Considerations
Temperature is an important consideration for making libraries:

* Keep libraries at temperatures ≤37°C, except where specifically noted.
* Place reagents on ice after thawing at room temperature.

### Equipment

* Review the programming instructions for your thermal cycler user guide to make sure that it is programmed appropriately using the heated lid function.
* It is acceptable to use the thermal cycler tracked heating lid function.


