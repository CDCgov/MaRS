June 8th, 2023; version 5.4.2; updated by [Marko Bajic](mailto:mbajic@cdc.gov)

> _Version 5.0: Major update made in library prep from NexteraXT to Flex protocol._

  * * *

**Table of contents**  

* [Introduction/General Overview](#intro)
* [Materials and Equipment](#materials)
* [Protocol Workflow](#workflow)
* [PET PCR Sample QC](#pet_pcr)
* [Gene PCR Enrichment](#gene_enrichment)
* [Electrophoresis](#electrophoresis)
* [PCR Amplicon Clean-Up](#pcr_cleanup)
* [Tagmentation of PCR Amplicon and Tagmentation Clean-Up​](#tagmentation)
* [Amplification of Tagmented DNA (Library Indexing)](#lib_indexing)
* [Library PCR Clean-Up](#lib_cleanup)
* [Library Pooling, Quantification, and Normalization](#lib_clustering)
* [Library Denaturing and MiSeq Sample Loading](#sample_loading)
* [Supporting Information](#supporting_info)


**IMPORTANT NOTICE**

This document provides information for an application for Illumina technology that has been demonstrated internally and may be of interest to external groups. This information is provided as‐is and is not an Illumina or CDC endorsed product and is not accompanied by any rights or warranties.

  * * *

<a id="intro"></a>

# Introduction #

**Standard Operating Procedure (SOP) describing how to prepare and sequence the full length _P. falciparum_ genes associated with antimalarial resistance on the Illumina MiSeq.**

* [_Plasmodium falciparum_ kelch 13 (_Pfk13_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_1343700)
* [_Plasmodium falciparum_ chloroquine resistant transporter (_Pfcrt_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0709000)
* [_Plasmodium falciparum_ multidrug resistant protein 1 (_Pfmdr1_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0523000)
* [_Plasmodium falciparum_ bifunctional dihydrofolate reductase-thymidylate synthase (_Pfdhfr_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0417200)
* [_Plasmodium falciparum_ dihydropteroate synthase (_Pfdhps_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0810800)
* [_Plasmodium falciparum_ cytochrome b (_Pfcytb_)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_MIT02300)

Human malaria is caused by six _Plasmodium_ species: *Plasmodium falciparum (Pf), P. vivax (Pv), P. malariae (Pm), P. ovale (Po) (P.o. curtisi and P.o. wallikeri)* and *P. knowlesi (Pk)*, which although zoonotic, can cause human infections in several South East Asian countries. Two of these, *Pf* and *Pv*, pose the greatest threat to global public health. About 3.2 billion people, half of the world's population, are at risk for malaria. In 2021, there were an estimated 247 million malaria cases in 84 malaria endemic countries, causing an estimated 619,000 deaths. In the U.S., an estimated 1,500 - 2,000 cases of malaria are imported annually. One of the greatest public health challenges for malaria control and elimination is the threat of drug resistant *Pf* parasites.

Previously effective anti-malarial treatments, chloroquine (CQ) and sulfadoxine/pyrimethamine (SP), are ineffective in many regions. Even more alarming, resistance to the least effective class of anti-malarial drugs, called artemisinins, has now emerged and spread in Southeast Asia, threatening malaria control and prevention programs globally.

Identifying and tracking drug resistance is critical for providing appropriate malaria prophylaxis and treatment guidelines. Molecular markers of resistance are available for several anti-malarial drugs, including artemisinins. Surveillance using molecular markers provides a robust system for the detection and tracking of resistant malaria parasites.

Below is a table of the major antimalarials and their associated with resistance molecular marker.

**Table 1. Antimalarials and associated resistance molecular markers**

|Antimalarial                            | Molecular Marker (Gene)  | Location     |
|----                                    |----                      |----          |
| Chloroquine                            | _Pfcrt_                    | Chr 7        |
| Artemisinin                            | _Pfk13_                    | Chr 13       |
| Atovaquone                             | _Pfcytb_                   | Mitochondria |
| Amodiaquine, lumefantrine, quinine     | _Pfmdr1_                   | Chr 5        |
| Pyrimethamine                          | _Pfdhfr_                   | Chr 4        |
| Sulfadoxine                            | _Pfdhps_                   | Chr 8        |
> Chr = chromosome; _Pfcrt_ = _Plasmodium falicparum_ chloroquine resistance transporter; _Pfk13_ = _Plasmodium falicparum_ kelch 13; _Pfcytb_ = _Plasmodium falicparum_ cytochrome b; _Pfmdr1_ = _Plasmodium falicparum_ multidrug resistance protein 1; _Pfdhfr_ = _Plasmodium falicparum_ bifunctional dihydrofolate reductase thymidylate synthase; _Pfdhps_ = _Plasmodium falicparum_ dydroxymethyldihydropterin pyrophosphokinase-dihydropteroate synthase;


  * * *

<!-- below page break can be used to define page breaks for pandoc conversion since it uses latex -->
<!-- \pagebreak -->

<br>

<a id="materials"></a>

#### Materials and Equipment ####

Please ensure all the necessary user‐supplied consumables and equipment are available before proceeding to sample preparation.

**Table 2. User‐Supplied Consumables**

| Consumable | Supplier |
| --- | --- |
| **Non-powdered sterile gloves** | General lab supplier |
| **Laboratory coat** | General lab supplier |
| **1.7 mL microcentrifuge tubes** | General lab supplier |
| **10 µL barrier pipette tips** | General lab supplier |
| **10 µL multichannel pipettes** | General lab supplier |
| **10 µL single channel pipettes** | General lab supplier |
| **20 µL barrier pipette tips** | General lab supplier |
| **20 µL multichannel pipettes** | General lab supplier |
| **20 µL single channel pipettes** | General lab supplier |
| **200 µL barrier pipette tips** | General lab supplier |
| **200 µL multichannel pipettes** | General lab supplier |
| **200 µL single channel pipettes** | General lab supplier |
| **1000 µL barrier pipette tips** | General lab supplier |
| **1000 µL multichannel pipettes** | General lab supplier |
| **1000 µL single channel pipettes** | General lab supplier |
| **PCR grade water** | General lab supplier |
| **RNase/DNase‐free 8‐well PCR strip tubes and caps** | General lab supplier |
| **[Optional] Disposable Polystyrene Reservoirs** | General lab supplier (Thomas Scientific Catalog #55501080) |
| **2X ABI TaqMan environmental buffer w/ Rox dye** | Applied Biosystems Catalog #4396838 |
| **Strip tubes 8X** | Agilent Catalog #410022 |
| **Strip Optical caps 8X** | Agilent Catalog #410024 |
| **Hardshell® 96‐well PCR Plates, clear semi skirted** | Bio-Rad Catalog #HSS9601 |
| **Phusion® High-Fidelity DNA Polymerase** | NEB Catalog #M0530L |
| **Deoxynucleotide (dNTP) Solution Mix** | NEB Catalog #N0447L |
| **Amplicon PCR Forward Primer (Standard desalting)** | (See tables 2-7) |
| **Amplicon PCR Reverse Primer (Standard desalting)** | (See tables 2-7) |
| **Lonza SeaKem® LE Agarose** | Lonza Catalog #50004 |
| **Nucleic Acid gel stain** | Biotum GelRed™ Nucleic Acid Gel Stain |
| **DNA gel loading dye** | Yakva Scientific 6X Orange-G Gel Loading Buffer #YSG |
| **Quick-Load 1 kb DNA ladder** | NEB Catalog #N0468 |
| **UltraPure™ 10X TBE Buffer** | Fisher Scientific Catalog #15581-044 |
| **AMPure XP beads for PCR Purification** | Beckman Coulter Life Sciences, Catalog #A63881 |
| **Illumina DNA Prep library kit** | Illumina, Catalog #20018705 (96 samples), or #20018704 (24 samples) |
| **\*IDT® for Illumina® DNA UD Index kits (plate)** | Illumina, Catalog #20027213 (Index set A), #20027214 (Index set B), #20027215 (Index set C), and #20027216 (Index set D).|
| **200 Proof Ethanol** | Decon Labs, Inc. Catalog #2716 |
| **AlumaSeal II aluminum seals** | Excel Scientific, Inc. Catalog #AF100 |
| **Clear, 8-strip PCR tubes domed caps** | LabSource, Catalog #T54-203-CS/10PK MFG# - 321-10-062 |
| **96‐well U-Shaped-Bottom Microplate** | Fisher Scientific, Catalog #7-200-720 |
| **Qubit® dsDNA HS Assay Kit** | Life Technologies Corporation Catalog #Q32854 |
| **Qubit™ Assay Tubes** | Thermo Fisher Scientific Catalog #Q32856 |
| **Agilent High Sensitivity D5000 ScreenTape** | Agilent Technologies, Catalog #5067-5592 |
| **Agilent High Sensitivity D5000 Reagents** | Agilent Technologies, Catalog #5067-5593 |
| **Agilent High Sensitivity D5000 Ladder** | Agilent Technologies, Catalog #5067-5594 |
| **1N NaOH** | Sigma-Aldrich, Inc. Catalog #SX0607H-6 |
| **Tris Hydrochloride, 1M Solution (pH 7.0/Mol. Biol.)** | Thermo Fisher Scientific Catalog #BP1756-100 |
| **MiSeq Reagent Kit V2 500 cycle kit** | Illumina Catalog #MS-102-2003 |
| **MiSeq Reagent V2 Nano Kit 500 cycle** | Illumina Catalog #MS-103-1003 |

**If you plan to pool >96 samples, you will need the Index Kit Set A, B, C, and D to provide unique multiplex combinations of indices**


**Table 3. User‐Supplied Equipment**

| Equipment | Supplier |
| --- | --- |
| **2-8°C Refrigerator** | General lab supplier |
| **-20°C Refrigerator** | General lab supplier |
| **Vortex** | General lab supplier |
| **4x Eppendorf PCR Cooler, iceless cold storage system for 96 well plates and PCR tubes** | Sigma-Aldrich, Inc (Z606634-1EA) |
| **Agilent ABI7500 or equivalent real-time PCR machine** | Agilent Technologies, Catalog #4351106 |
| **96‐well thermal cycler (with heated lid)** | General lab supplier |
| **Electrophoresis rig** | General lab supplier |
| **Magnetic stand‐96 or**  **96S Super Magnet Plate** | Life Technologies, Catalog #AM10027 or Alpaqua SKU A001322|
| **Microplate centrifuge** | General lab supplier |
| **Qubit 3.0 Fluorometer (or equivalent)** | Life Technologies Corporation, Catalog #Q33216 |
| **Agilent D4150 TapeStation System (or equivalent)** | Agilent Technologies, Catalog #G2992AA |
| **MiSeq Desktop Sequencer**| Illumina Inc. |

<a id="workflow"></a>


  * * *


#### Protocol Workflow ####

**NOTE:** The hands-on times are based on using 96-well format plates for each step.

1. [**PET-PCR Sample Quality Check**](#pet_pcr)<br />
Real-time PCR hands-on time 30 min / 96 samples; Cycle time 1.2 hours<br />
Reagents: Primers, 2X ABI TaqMan buffer, DNase PCR free water


2. [**PCR Reaction to Generate Amplicons**](#gene_enrichment)<br />
PCR hands-on time 30 min / 96 samples; Cycle time 2.5 hours<br />
Reagents: 10 µM Primers, HF Phusion Taq, 5X GC Buffer, 10 mM dNTPs, DNase PCR free water


3. [**Analysis of PCR Amplicons**](#electrophoresis)<br />
PCR amplicon electrophoresis hands-on time 10 min / 8 samples; Gel running time 30 min<br />
Reagents: Agarose, DNA loading dye, 1 kb DNA ladder, 1X TBE Buffer<br />
_If <20 samples, run all samples on the gel; If >20 samples, pick 20 samples with varying CT values and run on the gel_

4. [**PCR Amplicons Clean-Up**](#pcr_cleanup) <br />
Hands on time 30 min / 96 samples; Total time 40+ min / 96 samples<br />
Reagents: AMPure XP beads, fresh 70% EtOH, Nuclease-free water

5. [**Tagment Genomic DNA and Tagmentation Clean-Up**](#tagmentation) <br />
Hands on time 60 min / 96 samples; Total time 34 min / 8 samples<br />
Reagents: BLT, TB1, TSB, TWB<br />
_[optional]_ To assess tagmentation, run 1 µL sample on Agilent Bioanalyzer 2X and/or TapeStation 2X using High Sensitivity DNA chip

6. [**Amplification of Tagmented DNA (Library Indexing)**](#lib_indexing) <br>
Hands on time 35 min / 96 samples; Cycle time 38 min / 96 samples<br />
Reagents: EPM, Nuclease-free water, Index 1 and 2 primers

7. [**Library PCR Clean-Up**](#lib_cleanup) <br />
Hands on time 30 min / 96 samples; Total time 40+ min / 96 samples<br />
Reagents: SPB, RSB, Nuclease-free water, fresh 80% EtOH

8. [**Library Pooling, Quantification, and Normalization**](#lib_clustering) <br />
Hands on time 30+ min / 96 samples; Total time 40+ min / 96 samples<br />
Reagents: Sample Buffer, D5000 Ladder, ScreenTape; Qubit dsDNA HS Buffer and Reagent, Standard #1 and #2

9. [**Library Denaturing and MiSeq Sample Loading**](#sample_loading) <br />
Hands on time 30 min / pooled samples; Total time 30 min / pooled samples<br />
Reagents: Resuspension Buffer, HT1, 0.2N NaOH, PhiX Control Kit v3, 200 mM Tris-HCl pH7.0

10. **Analysis of NGS data** <br />
Hands on time 5 min / 96 samples; Total time 15-25 min / 96 samples<br />
Method: [MaRS analysis pipeline](github.com/CDCgov/MaRS/tree/master/Geneious_workflow)

11. **Standardized SNPs Reports Generated**



  * * *


<a id="pet_pcr"></a>

## Sample QC ##

This step uses a real time PCR assay, [PET-PCR](https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0179178),  to assess the quality and quantity of starting DNA material. The readout includes an estimation of _all_ DNA in the sample, host and parasite.

- **NOTE:** This assay is used to identify which samples should be included and/or excluded for downstream procedures. Any sample with a CT value greater than 34, should be excluded. We recommend grouping all samples with a CT > 34 and performing separate amplicon PCRs and electrophoresis for these samples. If the electrophoresis yields positive results (e.g. positive and clear bands on gel) you may procced with downstream procedures. Samples with a CT > 34 have lower  parasite DNA and subsequent gene amplification success rates are variable.

**Consumables**

**Table 4. PET-PCR Consumables**

| Item | Quantity | Storage |
| --- | --- | --- |
| **Primers – FAM labled genus and HEX labeled _falciparum_ (see below)** | 0.25-0.5 µL per sample | 2° to 8°C |
| **TaqMan 2X Environmental buffer** | 10 µL per sample | 2° to 8°C |
| **Nuclease-free water** | 6.25 µL per sample | Room temperature |
| **Strip tubes 8X** | Up to 8 samples per strip | Room temperature |
| **Strip Optical caps 8X** | Up to 8 samples per strip | Room temperature |

**Preparation**

- All stock primers should be prepared at a 10 µM concentration.<br>
- DNA samples should be stored at 4&deg;C until testing or -20&deg;C for long term storage.<br>
- Store all primer stocks at -20&deg;C for up to 1 year.<br>
- Unopened tubes of ABI TaqMan Environmental Buffer should be stored at -20&deg;C for a maximum of six months. Once thawed, store at 4&deg;C for up to six months. The reagent must be used within the expiration date provided by the manufacturer.<br>
- All samples should be tested in duplicates or triplicates in some special cases (e.g. very low density situations).<br>
- Typically, the genus _P. falciparum_ multiplex assay should be run first on all samples. All genus-positive sample are subsequently tested for _P. ovale, P. malaria_ and _P. vivax_, in order to determine the species.

**Procedure**

**Initial Set up**

- The PET-PCR reaction mix is prepared by mixing the TaqMan environmental buffer, primers, and water as shown below.
- Determine the number of reactions you need to run by multiplying the total number of samples you have to test (including your positive and negative controls) by two because every sample will be tested in duplicates.
  - For example, if you are testing 10 samples, you will multiple this by two to give you 20. Add two extra reactions to account for loss of solution during pipetting. This gives you a total number of 22 reactions. Multiply this number with the volumes below for each component to give you the total master-mix volume required for your experiment.
- In a 1.5 mL tube, prepare your master-mix by multiplying the volumes shown below with the total number of reactions you need to run (e.g. 22 as described above).


**Primers and PCR Conditions**

The table below shows the Genus and _P. falciparum_ primers and PCR conditions for a multiplex reaction:  

**Table 5. Multiplexing Genus and _P. falciparum_ species specific primers at 10 µM**

| Master mix              | Reaction volume | x samples + 10% | Final [conc] |
|----                     |----             |----             |----          |
|  Water                  | 6.25 µL         |                 |              |
|  2X ABI TaqMan buffer   | 10.00 µL        |                 |    1x        |
|  Genus F                | 0.50 µL         |                 |    0.250 µM  |
|  FAM-genus R            | 0.50 µL         |                 |    0.250 µM  |
|  P.f. F                 | 0.50 µL         |                 |    0.250 µM  |
|  HEX-P.f. R             | 0.25 µL         |                 |    0.125 µM  |
|  **TOTAL**              | **18.0 µL**     |                 |              |
|  **Add last**           |                 |                 |              |
|  DNA                    | 2.0 µL          |                 |              |

>Primer (5' to 3'):<br />
>**Genus 18sFor (5' to 3'):** 5'-GGCCTAACATGGCTATGACG-3'<br />
>**Genus 18sRev (5' to 3'):** 5'-aggcgcatagcgcctggCTGCCTTCCTTAGATGTGGTAGCT-3' (**FAM-labeled: based on the 18s rRNA gene)**<br />
>_**P. falciparum**_ **For (5' to 3'):** 5'-ACCCCTCGCCTGGTGTTTTT-3'<br />
>_**P. falciparum**_ **Rev (5' to 3'):** 5'-aggcggataccgcctggTCGGGCCCCAAAAATAGGAA-3' **(HEX-labeled: based on the r364 target)**<br />

**Thermocyclying conditions**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 95&deg;C                           | 15:00      |  
|        |                                    |            |
| **2**  | **95&deg;C**                       | **0:20**   |      
| **3**  | **63&deg;C**                       | **0:40**   |      
| **4**  | **72&deg;C**                       | **0:30**   |    
|        | **Repeat Steps 2-4 for 44 cycles<br /> (45 total)** |            |  
|        |                                    |            |
| 5      | 4&deg;C                            | Infinity   |

**Adding the DNA Samples**

1. Mix the prepared master-mix well by vortexing briefly.

2. Centrifuge the tubes for 5 seconds to remove any solution trapped in the cap.

3. Arrange the optically clear PCR tubes on a PCR-tube rack following the PCR sample sheet. Add 18 µL of the PET-PCR master mix prepared above to each PCR well. Loosely put on the lids of the wells filled with master mix solution.

4. Return all reagents to the freezer and refrigerator before proceeding to the next step.

5. Take the assembled plate containing the tubes with PCR master mix solution to the PCR template area.

6. Add 2 µL of the unknown DNA samples to the wells with the master-mix according to the sample sheet. Cap the well tightly after adding the sample. The total volume of PCR reaction is 20.0 µL after addition of the template.

7. Add positive control DNA to each positive control well with master-mix. Cap the wells after each positive control is added.

8. Add 2.0 µL of DNase-free water to the wells designated as the no-template control (NTC) and close that well tightly.

9. Make sure each sample has been added to the correct well and that all wells are tightly capped.

10. Briefly centrifuge your strip tubes to remove any solution trapped on the walls of the wells.

11. Make sure there are no bubbles in the well.

- **NOTE:** _The recommended minimum amount of template DNA is 2.0 µL. This can be adjusted appropriately depending on the sample parasitemia._

**PCR-Cycling Parameters**

1. Start the real-time PCR thermocycler according to the manufacturer's guidelines.

2. Program the software to detect fluorescence through FAM, HEX and ROX filters all wells. ROX is to be detected as a reference dye.

3. Program the software to run the cycling conditions shown under table 5.

4. Fluorescence data should be collected at the amplification plateau.

**Interpreting Results**

1. Use default settings to set the dye thresholds.  

2. If the calculated thresholds are located within the background noise, they should be manually set to a level slightly higher than the background. Such alterations should be done with only one dye displayed at the time.

3. Positive specimens are those that yield a fluorescence signal above the threshold value in the wells where samples or controls were loaded
  * Positive PCR: A positive sample produces a fluorescence signal above the threshold/noise level. Positive samples are designated a Ct value below 40.0.
  * Negative PCR: No fluorescence signal above the threshold/noise level. Negative samples have no Ct or have a Ct value above 40.0.

- **NOTE:** _The negative controls must be negative (no Ct or above 40.0). The positive controls must be positive (designated by Ct value below 40.0). The test should be repeated if the NTC has a positive Ct value, or if the positive control yields no positive results._

**For more information, please see:**
[Lucchi, N.W., et al., _Molecular diagnosis of malaria by photo-induced electron transfer fluorogenic primers: PET-PCR._ PLoS One, 2013. 8 (2): p. e56677.](https://pubmed.ncbi.nlm.nih.gov/23437209/)


  * * *


## Gene Enrichment & QC ##

#### [I: Gene PCR Amplification](#gene_amp) ####
#### [II: QC by Electrophoresis](#electrophoresis) ####

### Recommended Positive Controls and Expected SNPs ###

We highly recommend using at least **three positive controls** with known SNP profiles. These will be first analyzed to confirm the known SNPs and ensure the rest of the sequencing run was successful.

We routinely use the following controls:

| Control strain | _Pfcrt_         | _Pfmdr1_        | _Pfdhfr_        | _Pfdhps_    | _Pfk13_       |
| ---            | ---           | ---           | ---           | ---       | ---         |
| 7G8            | **S**VMN**T** | NED**FCDY**   | C**I**C**N**I | S**G**KAA | _wild type_ |
| DD2            | CV**IET**     | **Y/F**EDFCDY | C**IRN**I     | S**G**KAA | _wild type_ |
| HB3            | CVMNK         | **Y/F**EDFCDY | C**I**C**T**I | S**G**KAA | _wild type_ |

> **BOLD** indicates codon position with mutations<br />
> Codon positions:_Pfcrt_ :72-76; _Pfmdr1_ :86,130, 144, 184, 1034, 1042, 1109, 1246; _Pfdhfr_ : 50, 51, 59, 108, 164; _Pfk13_ : 18 - 715

Controls can be ordered directly through BEI: https://www.beiresources.org/


<a id="gene_enrichment"></a>

#### I: Gene PCR Amplification ####

This step uses PCR to amplify template from a DNA sample using region of interest-specific primers.

User‐defined forward and reverse primers are used to amplify templates from genomic DNA. A subsequent limited‐cycle amplification step is performed to add multiplexing indices and Illumina sequencing adapters. Libraries are normalized and pooled, and sequenced on the MiSeq system using v2 reagents.

**Preperation**

**Initial Set up**

- Ensure that the No-DNA and DNA-only UV stations have all the appropriate pipettes and tip sizes.
- Clean up all pipettes and lab bench area using using 10% bleach followed by 70% ethanol.
- Turn on the No-DNA (PCR master mix) and DNA-only UV station for 30 minutes.
- If you have not already done so, create a 10 µM working stock solution of your primers, using C1V1 = C2V2 to calculate the appropriate volume needed to make a working stock solution.
- Label all freshly made and newly opened items with the date and your name initials.
- Get the appropriate number of PCR plates and/or PCR tubes and place them in the no-DNA UV station.
- Label the PCR tubes (use a printout template for PCR plates) with sample IDs.
- Get the appropriate number EppendorfPCR Cooler plates from the freezer, wipe down with 70% ethanol and place them in the UV hoods (one or more sets each in the no-DNA and DNA-only UV stations). Turn on the UV stations for another 30 min.
- Be sure to reserve the appropriate number of thermocyclers and have the appropriate cycling conditions set up.

**Pre-PCR master mix prep**

- Let Primers, dNTPs, and GC Buffer defrost at room temperature (10-15 min). Once defrosted, mix gently _(vortexed)_ and centrifuge briefly prior to use.
- **_DO NOT thaw and/or vortex or mix the HF Taq._**
- All PCR reactions **must be assembled on the** Eppendorf PCR Cooler plates
- **Always add the Taq last when making your master mix and DO NOT vortex and/or pipette after adding Taq.**
- If you forget to return any of the reagents, especially the Taq, to its appropriate storage conditions (i.e., leave it out at room temperature), record the date and time of when it happened in your lab notebook, and discard the Taq.

**PCR plate set up**
- Calculate appropriate volumes for master mix based on number of samples to be included in reaction; multiply each reagent volume times the total number of samples + 10% total number of samples (to have enough mix for all samples)
- Final volume of master mix is given in **Tables 6.1 - 6.6**
- **NOTE:**  _If the number of samples is >5, make a master mix for at least 6 samples to avoid not having enough mix._

**Procedure**

1. Set up PCR reaction consisting of water, GC Buffer, dNTPs, primers, [Phusion High-Fidelity DNA Polymerase](https://www.neb.com/products/m0530-phusion-high-fidelity-dna-polymerase#Product%20Information), and DNA in the order given in **Tables 6.1 – 6.6:**

2. Seal plates and/or PCR tubes.

3. Once tubes and/or plates are sealed, keep them in the Eppendorf PCR Cooler plates.

4. **Pre-heat** the thermal cycler to 98&deg;C prior to placing PCR plates and/or PCR tubes into the thermal cycler. Pre-heating to 98&deg;C should take 0:30 of the 3:00 min.


**Primers and PCR Conditions**

The tables below show primers and PCR conditions for the following antimalarial drug resistance associated genes: _Pfcrt_ [6.1](#6.1), _Pfk13_ [6.2](#6.2), _mitochondria_ [6.3](#6.3), _Pfcytb_ [6.3a](#6.3a), _Pfdhps_ [6.4](#6.4), _Pfdhfr_ [6.5](#6.5), _Pfmdr1_ [6.6](#6.6). Genes for parasite fingerprinting: _Pfs47_ [6.7](#6.7), and _Pfcpmp_ [6.8](#6.8)

**IMPORTANT**: While the master mix conditions will be the same for all genes, the thermocycling conditions will differ, specifically the annealing temperatures.  

<a id="6.1"></a>
**Table 6.1. Gene: _Pfcrt_ (3,109 bp); <ins>Primers at 10 µM</ins>**

| Master mix      | Reaction volume | x samples + 10% | Final [conc] |
|----             |----             |----             |----          |
|  5X GC Buffer   | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)  | 1.0 µL          |                 |    0.2 nM    |
|  _mars_crt-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_crt-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water          | 32.0 µL         |                 |              |
|  **Add last**   |                 |                 |              |
|  HF Phusion Taq | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**      | **46.0 µL**     |                 |              |
|  **Add**        |                 |                 |              |
|  Template DNA   | 4.0 µL per well |                 |              |
|  **TOTAL**      | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_crt-fwd_: _TTACATATAACAAAATGAAATTCGC_<br />
>_mars_crt-rev_: _TATTGTGTAATAATTGAATCGACG_

**Thermocyclying conditions for _Pfcrt_; <ins> Primers at 10 µM</ins>**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **62&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |

<a id="6.2"></a>
**Table 6.2.  _Pfk13_ (2,181 bp); <ins> Primers at 10 µM</ins>**

| Master mix      | Reaction volume | x samples + 10% | Final [conc] |
|----             |----             |----             |----          |
|  5X GC Buffer   | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)  | 1.0 µL          |                 |    0.2 nM    |
|  _mars_k13-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_k13-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water          | 32.0 µL         |                 |              |
|  **Add last**   |                 |                 |              |
|  HF Phusion Taq | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**      | **46.0 µL**     |                 |              |
|  **Add**        |                 |                 |              |
|  Template DNA   | 4.0 µL per well |                 |              |
|  **TOTAL**      | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_k13-fwd_: _CTATGACGTATGATAGGGAATCTGG_<br />
>_mars_k13-rev_: _CTGGGAACTAATAAAGATGGGCC_

**Thermocyclying conditions for _Pfk13_; <ins> Primers at 10 µM</ins>**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **58&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


<a id="6.3"></a>
**Table 6.3. _Mitochondria_ (5,967 bp); <ins> Primers at 10 µM</ins>**

**NOTE:** If experiencing issues with amplifying the full-length mitochondrial genome, consider amplifying only the _cyt-b_ gene instead for characterizing molecular markers associated with Malarone (atovaquone/proguanil) resistance. See Table 6.3a below.

| Master mix      | Reaction volume | x samples + 10% | Final [conc] |
|----             |----             |----             |----          |
|  5X GC Buffer   | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)  | 1.0 µL          |                 |    0.2 nM    |
|  _mars_mit-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_mit-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water          | 32.0 µL         |                 |              |
|  **Add last**   |                 |                 |              |
|  HF Phusion Taq | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**      | **46.0 µL**     |                 |              |
|  **Add**        |                 |                 |              |
|  Template DNA   | 4.0 µL per well |                 |              |
|  **TOTAL**      | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_mit-fwd_: _AAGCTTTTGGTATCTCGTAAT_<br />
>_mars_mit-rev_: _TATTATAATATAACTCTACAAAGTTGAAC_

**Thermocyclying conditions for _Mitochondria_; <ins> Primers at 10 µM</ins>**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **50&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **6:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


<a id="6.3a"></a>
**Table 6.3a. _Cytochrome b_ (937 bp); <ins> Primers at 10 µM</ins>**

| Master mix       | Reaction volume | x samples + 10% | Final [conc] |
|----              |----             |----             |----          |
|  5X GC Buffer    | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)   | 1.0 µL          |                 |    0.2 nM    |
|  _mars_cytb-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_cytb-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water           | 32.0 µL         |                 |              |
|  **Add last**    |                 |                 |              |
|  HF Phusion Taq  | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**       | **46.0 µL**     |                 |              |
|  **Add**         |                 |                 |              |
|  Template DNA    | 4.0 µL per well |                 |              |
|  **TOTAL**       | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_cytb-fwd_: _CTATTAATTTAGTTAAAGCACAC_<br />
>_mars_cytb-rev_: _ACAGAATAATCTCTAGCACCA_

**Thermocyclying conditions for _cytochrome b_; <ins> Primers at 10 µM</ins>**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **60&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **3:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


<a id="6.4"></a>
**Table 6.4. _Pfmdr1_ (4,155 bp); <ins> Primers at 10 µM</ins>**

| Master mix       | Reaction volume | x samples + 10% | Final [conc] |
|----              |----             |----             |----          |
|  5X GC Buffer    | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)   | 1.0 µL          |                 |    0.2 nM    |
|  _mars_mdr1-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_mdr1-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water           | 32.0 µL         |                 |              |
|  **Add last**    |                 |                 |              |
|  HF Phusion Taq  | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**       | **46.0 µL**     |                 |              |
|  **Add**         |                 |                 |              |
|  Template DNA    | 4.0 µL per well |                 |              |
|  **TOTAL**       | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_mdr1-fwd_: _TGGTAACCTCAGTATCAAAG_<br />
>_mars_mdr1-rev_: _CATCTTGTGCTGATAATAATTC_

**Thermocyclying conditions for _Pfmdr1_; <ins> Primers at 10 µM**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **60&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |  
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


<a id="6.5"></a>
**Table 6.5. _Pfdhfr_ (2,067 bp); <ins> Primers at 10 µM</ins>**

| Master mix       | Reaction volume | x samples + 10% | Final [conc] |
|----              |----             |----             |----          |
|  5X GC Buffer    | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)   | 1.0 µL          |                 |    0.2 nM    |
|  _mars_dhfr-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_dhfr-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water           | 32.0 µL         |                 |              |
|  **Add last**    |                 |                 |              |
|  HF Phusion Taq  | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**       | **46.0 µL**     |                 |              |
|  **Add**         |                 |                 |              |
|  Template DNA    | 4.0 µL per well |                 |              |
|  **TOTAL**       | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_dhfr-fwd_: _TTTTTACTAGCCATTTTTGTATTCC_<br />
>_mars_dhfr-rev_: _TTAACCGTTCAGGTAATTTTGTCA_

**Thermocyclying conditions for _Pfdhfr_; <ins> Primers at 10 µM**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **58&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


>Primers adapted from: [SC, Carlton JM. 2016. A Method for Amplicon Deep Sequencing of Drug Resistance Genes in Plasmodium falciparum Clinical Isolates from India. J Clin Microbiol 54:1500–1511.](https://pubmed.ncbi.nlm.nih.gov/27008882/)


<a id="6.6"></a>
**Table 6.6. _Pfdhps_ (2,817 bp); <ins> Primers at 10 µM</ins>**

| Master mix       | Reaction volume | x samples + 10% | Final [conc] |
|----              |----             |----             |----          |
|  5X GC Buffer    | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)   | 1.0 µL          |                 |    0.2 nM    |
|  _mars_dhps-fwd_ | 1.25 µL         |                 |    0.25 µM   |
|  _mars_dhps-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water           | 32.0 µL         |                 |              |
|  **Add last**    |                 |                 |              |
|  HF Phusion Taq  | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**       | **46.0 µL**     |                 |              |
|  **Add**         |                 |                 |              |
|  Template DNA    | 4.0 µL per well |                 |              |
|  **TOTAL**       | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'): <br />
>_mars_dhps-fwd_: _AATATTTGCGCCAAACTTTTTA_<br />
>_mars_dhps-rev_: _TTTATTTCGTAATAGTCCACTTTTGAT_

**Thermocyclying conditions for _Pfdhps_; <ins> Primers at 10 µM**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            |
| **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **58&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


>Primers adapted from: [SC, Carlton JM. 2016. A Method for Amplicon Deep Sequencing of Drug Resistance Genes in Plasmodium falciparum Clinical Isolates from India. J Clin Microbiol 54:1500–1511.](https://pubmed.ncbi.nlm.nih.gov/27008882/)


<a id="6.7"></a>
**Table 6.7. _Pfs47_ (1,320 bp); <ins> Primers at 10 µM</ins>**

| Master mix       | Reaction volume | x samples + 10% | Final [conc] |
|----              |----             |----             |----          |
|  5X GC Buffer    | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)   | 1.0 µL          |                 |    0.2 nM    |
| _mars_pfs47-fwd_ | 1.25 µL         |                 |    0.25 µM   |
| _mars_pfs47-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water           | 32.0 µL         |                 |              |
|  **Add last**    |                 |                 |              |
|  HF Phusion Taq  | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**       | **46.0 µL**     |                 |              |
|  **Add**         |                 |                 |              |
|  Template DNA    | 4.0 µL per well |                 |              |
|  **TOTAL**       | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_pfs47-fwd_: _ATGTGTATGGGAAGAATGATCAG_<br />
>_mars_pfs47-rev_: _TCATATGCTAACATACATGTAAAAAATTAC_

**Thermocyclying conditions for _Pfs47_; <ins> Primers at 10 µM**:

| Step   | Temperature                        | Time (min) |
| :--:   |:----:                              | :----:     |
| 1      | 98&deg;C                           | 3:00       |  
|        |                                    |            | **2**  | **98&deg;C**                       | **0:30**   |      
| **3**  | **58&deg;C**                       | **0:30**   |      
| **4**  | **65&deg;C**                       | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |  
|        |                                    |            |
| 5      | 65&deg;C                           | 10:00      |
| 6      | 4&deg;C                            | Infinity   |


<a id="6.8"></a>
**Table 6.8. _Pfcpmp_ (425 bp); <ins> Primers at 10 µM</ins>**

| Master mix      | Reaction volume | x samples + 10% | Final [conc] |
|----             |----             |----             |----          |
|  5X GC Buffer   | 10.0 µL         |                 |    1x        |
|  dNTPs (10 mM)  | 1.0 µL          |                 |    0.2 nM    |
| _mars_cpmp-fwd_ | 1.25 µL         |                 |    0.25 µM   |
| _mars_cpmp-rev_ | 1.25 µL         |                 |    0.25 µM   |
|  Water          | 32.0 µL         |                 |              |
|  **Add last**   |                 |                 |              |
|  HF Phusion Taq | 0.5 µL          |                 |    1 unit    |
|  **TOTAL**      | **46.0 µL**     |                 |              |
|  **Add**        |                 |                 |              |
|  Template DNA   | 4.0 µL per well |                 |              |
|  **TOTAL**      | **50.0 µL**     |                 | 50 - 250 ng  |

>Primers (5'to 3'):<br />
>_mars_cpmp-fwd_: _GTCATTAAAATTTATGGATTATATATGG_<br />
>_mars_cpmp-rev_v2:_GTTACTATCAAGATCGTTAATATC_

**Thermocyclying conditions for _Pfcpmp_; <ins> Primers at 10 µM**:

| Step   | Temperature                                   | Time (min) |
| :--:   |:----:                                         | :----:     |
| 1      | 98&deg;C                                      | 3:00       |  
|        |                                               |            |
| **2**  | **98&deg;C**                                  | **0:30**   |      
| **3**  | **54&deg;C**                                  | **0:30**   |      
| **4**  | **65&deg;C**                                  | **5:00**   |    
|        | **Repeat Steps 2-4 for 29 cycles <br />(30 total)** |            |  
|        |                                               |            |
| 5      | 65&deg;C                                      | 10:00      |
| 6      | 4&deg;C                                       | Infinity   |


**SAFE STOPPING POINT**
If you do not immediately proceed to Electrophoresis, seal plate with adhesive seal and store it at 2° to 8°C for up to a week.

 * * *

<a id="electrophoresis"></a>

#### II. QC by Electrophoresis ####

This step is necessary to ensure successful amplification of amplicons. It is recommended to run at least 25% of the total samples, all no-template and negative controls on the gel to confirm amplification was successful and no contamination occurred. Please note PCR amplification can be affected by numerous factors, including but not limited to, DNA quality and quantity.

- **First run only the positive, no-template and negative controls**.
- If unexpected results are obtained, record the time and information in your notebook, then toss the samples. **Amplification will need to be repeated** [I: Gene PCR Amplification](#gene_amp)

- If expected results are obtained, move on to running at minimum 25% of the samples along with the positive, NTC and NC.

Of the 25% total samples, ensure to select representative samples of varying parasitemia levels (based on PET-PCR CTs or microsopy data).

**Consumables**

**Table 7a. Electrophoresis Consumables**

| Item | Quantity | Storage |
| --- | --- | --- |
| **Agarose** | 1g (for a 1% gel) | Room temperature |
| **1x solution of 10X TBE Buffer and deionized water** | 100 mL (for a 1% gel) | Room temperature |
| **Nucleic Acid Gel Stain** | 5 µL per 100 mL of buffer | Room temperature |
| **Orange Dye** | 2 µL per 8 µL PCR product | Room temperature |

**Preparation**

 * Create a 1X solution of 10X TBE Buffer and deionized water before you start.
 * For a 1% agarose gel, 1.0 gram of agarose + 100 mL of buffer will fill a medium gel chamber
 * For a 1.5% agarose gel, 1.5 gram of agarose + 100 mL of buffer will fill a medium gel chamber

**Procedure**

1. Choose an Erlenmeyer flask that is 2-4 times the volume of the solution and place a stirring rod into the flask.

2. Weigh the agar to the desired concentration.

3. Add the appropriate amount of buffer for the desired concentration.

4. Dissolve the agar in the microwave by heating the solution on high power until it comes to a boil. Watch the solution closely; **DO NOT** allow solution to boil over.

5. Remove the flask with pot holders and gently swirl to re-suspend any settled agar.

6. Repeat steps 4-5 until all the agar is dissolved (no transparent agarose clumps should be present).

7. Allow the solution to cool on a stirring plate until you can comfortably hold the flask with your hands.

8. Using a 10 µL pipette, add nucleic acid gel stain to the solution. For every 100 mL of buffer, add 5 µL of gel stain. Swirl solution to mix, making sure as little bubbles as possible are created.

9. Pour the cooled solution into the gel device and ensure no bubbles are present. Place the comb into the gel and allow the gel to sit undisturbed for at least 15 minutes or until the gel has become firm (the color will change from clear to slightly milky in color).

10. When gel has solidified, ensure the wells are aligned with the black (negative) nodes on the electrophoresis chamber and fill with buffer until it covers about a centimeter above the gel. Remove the comb by pulling it gently into an upward direction.

11. Combine a mixture of 2 µL of orange dye and 8 µL of each sample and load 8 µL of that mixture into each well.

12. Be sure to include a 1 kb and 100 bp reference dye ladders, one on each side of the gel (no orange dye necessary).

13. Place the lid on the chamber box and connect the black node to the negative terminal and the red node to the positive terminal. Turn on the power supply and adjust the voltage to 100-130 volts.

14. Run gel for 40-45 minutes; check the gel at 40 min; the samples should nearly reach the end of the gel. **DO NOT** allow samples to run off the gel.

15. Turn off the power supply, disconnect the electrodes, and remove the lid from the gel device.

16. Remove the gel from the chamber and take to the gel reading station for analysis.

17. Once amplification is confirmed, the samples can proceed to [PCR amplicon Clean-Up](#pcr_cleanup).


 * * *

<a id="pcr_cleanup"></a>

## PCR Amplicon Clean-Up ##

This step uses AMPure XP beads to clean up your PCR amplicon gene product(s). You can locate Agencourt AMPure XP PCR Purification Instructions for Use. **PLEASE SEE VENDOR PROTOCOL [HERE](https://www.beckman.com/search#q=A63881&t=coveo-tab-techdocs).**

**Consumables**

**Table 8. PCR Amplicon Purification Consumables**

| Item                                | **Quantity**               | Storage           |
| ---                                 | ---                        | ---               |
| AMPure XP beads                     | 90 µL per 50 µL of sample  | 2&deg;C to 8&deg;C|
| Freshly Prepared 70% Ethanol (EtOH) | 400 µL per sample          | Room temperature  |
| Nuclease free water; # 25-055-CM    | 40 µL per sample           | Room temperature  |
| 96‐well 0.2 mL PCR plate            | 1 plate                    |                   |
| [Optional] Microseal 'B' film       |                            |                   |
| 96‐well U-Shaped-Bottom Microplate  | 1 plate                    |                   |


**Preparation**

- Determine whether or not a plate transfer is necessary.  If the PCR reaction volume multiplied by 2.8 exceeds the volume of the PCR plate, a transfer to a 300 µL round bottom plate or a 1.2 mL deep-well plate is required.

- Bring the AMPure XP beads to room temperature- **wait at least 30 minutes**.

- Once at room temperature, shake the Agencourt AMPure XP bottle to re-suspend any magnetic particles that may have settled. Ensure magnetic beads are well (evenly) distributed before adding them to samples.

- **_Thaw RSB at room temperature._**

- **NOTE**: Proceed very cautiously during this procedure and take your time.

- **NOTE**: 70% ethanol is hygroscopic. That is, when opened the ethanol will both evaporate and absorb water over time. Re-use eventually will be at a lower concentration. There is also miscibility involved with ethanol and water. For example, measuring out 80 mL of ethanol and topping off to 100 mL with water will generate ~65% ethanol. Measuring 80 mL ethanol and 20 mL water separately, then combining them will generate ~95 mL of 70% ethanol. Make sure to use molecular biology grade water (DNAse, RNase and Protease free).


**Procedure**

1. Centrifuge the Library Amplification plates at 1,000 × g at 20°C for 1 minute to collect condensation, carefully remove seal.

2. Combine (pool) each PCR gene amplicon into a new 96‐well U-Shaped-Bottom Microplate using **table 7b** below.

    **Note:** The PCR efficiency for each of the genes varies (_Pfk13_ > _Pfdhfr_ > _Pfdhps_ > _Pfmdr1_ > _Pfcrt_; highest to lowest PCR efficiency). Thus, its important to add approximately the same total concentration of each PCR gene amplicon to each pool of gene amplicons. Using the gel from [QC by Electrophoresis](#electrophoresis), use the [gel analyzer](http://www.gelanalyzer.com/?i=1) program to determine estimated concentration of your genes. Based on this analysis, adjust the total volume added of each gene, making sure to *always add at minimum 5 µL* from each gene PCR amplicon to the final pool.

**Table 7b. Amounts for pooling PCR amplicons**

| Item            | **Quantity**         | Storage    |
| ---             | ---                  | ---        |
| 2 genes         | combine 25 µL each   | 2° to 8°C  |
| 3 genes         | combine 16 µL each   | 2° to 8°C  |
| 4 genes         | combine 12.5 µL each | 2° to 8°C  |
| 5 or more genes | combine 10 µL each   | 2° to 8°C  |

> This should yield a total of 50 µL of combined PCR gene product.


3. Gently shake the AMPure XP beads for 30 seconds to make sure that the beads are evenly dispersed. Add an appropriate volume of beads to a trough depending on the number of samples being processed and desired fragment selection. **The Illumina DNA Prep library kits typically yield <ins>insert sizes around the 500 bp range.</ins>**

   **NOTE** To maximize recovery of smaller fragments from the bead cleanup step, use the following conditions (AMPure XP volume calculations in the table below are based on 50 µL sample volume):

| Input Size (bp)    | AMPure XP Recommendation      | AMPure XP Volume (µL)  |
| ---                | ---                           | ---                    |
| < 300              | 1.8x AMPure XP                | 90.0                   |
| 300 - 500          | 1.8x AMPure XP                | 90.0                   |
| <ins>**500**</ins> | <ins>**0.6x AMPure XP**</ins> | <ins>**30.0**</ins>    |
| gDNA               | 0.6x AMPure XP                | 30.0                   |

4. Using a multichannel pipette, add an appropriate volume of beads per sample based on your input size. While adding the beads to the samples, gently pipette the entire volume up and down 10 times. Change tips between columns.

5. Incubate the mixed samples at room temperature for 5 minutes.

6. Place the library amplification plate on a magnetic stand for 2 minutes. WAIT for the solution to clear before proceeding.

7. With the library amplification plate on the magnetic stand, use a multichannel pipette to _carefully_ remove and discard all the supernatant. Change tips between samples.

    **DO NOT disturb the magnetic beads.**

8. With the library amplification plate on the magnetic stand, wash the beads with freshly prepared 70% ethanol as follows:

    - A. Using a multichannel pipette, add 200 µL of freshly prepared 70% ethanol to each sample well.

    - B. Incubate the plate on the magnetic stand for 30 seconds at room temperature.

    - C. Carefully remove excess ethanol using a P20 multichannel pipette.

      **NOTE: The beads are not drawn out easily when in alcohol, so it is not necessary to leave any supernatant behind.**

9. With the library amplification plate still on the magnetic stand, perform a second ethanol was as follows:

    - A. Using a multichannel pipette, add 180 µL of freshly prepared 70% ethanol to each sample well.

    - B. Incubate the plate on the magnetic stand for 30 seconds at room temperature.

    - C. Carefully remove excess ethanol using a P20 multichannel pipette.

    - D. Use a P20 multichannel pipette with fine pipette tips to remove excess ethanol.

10. With the library amplification plate still on the magnetic stand, allow the beads to air-dry for at least 3 minutes, then check every 2 minutes until no ethanol remains.

    **NOTE: Make sure not to over dry the beads. Bead pellets will appear cracked if over dried.**

11. Remove the library amplification plate from the magnetic stand. Using a multichannel pipette, add 32 µL nuclease free water (# 25-055-CM, Cell Culture Grade Water, Sterile) to each well of the library amplification plate. While adding the nuclease free water to the beads, gently pippete mix the beads and water up and down 10 times. Change tips after each column.

12. Incubate at room temperature for 2 minutes.

13. Place the plate back on the magnetic stand for 2 minutes or until the supernatant has cleared.

14. Using a multichannel pipette, carefully transfer 30 µL of the supernatant from the Library amplification plate to a new 96‐well PCR plate. Change tips between samples to avoid cross‐contamination.


**SAFE STOPPING POINT**

If you do not plan to proceed to *Tagment Genomic DNA*, seal the plate with Microseal "B" adhesive seal. Store the plate at ‐15° to ‐25°C for up to a week.

  * * *

## NGS Library Prep ##

#### [I: Tagment Genomic DNA](#tagmentation) ####
#### [II: Post Tagmentation Clean-Up](#tag_cleanup) ####
#### [III: Amplification of Tagmented DNA (Library Indexing)](#lib_indexing) ####
#### [IV: NGS Library Clean‐Up](#lib_cleanup") ####

<a id="tagmentation"></a>

#### I: Tagment Genomic DNA ####

The tagmentation step uses the Bead-Linked Transposomes (BLT) to tagment DNA. This process fragments and tags the DNA with adapter sequences. The Post Tagmentation Clean up step washes the adapter-tagged DNA on the BLT before PCR amplification. See [How on-bead tagementation works](https://www.illumina.com/techniques/sequencing/ngs-library-prep/tagmentation.html)


**Consumables**

**Table 9. Tagment Genomic DNA Consumables**

| Item                               | Quantity           | Storage                |
| ---                                | ---                | ---                    |
| **TB1 (Tagmentation Buffer 1)**    | 10 µL per sample   | -15&deg;C to -25&deg;C |
| **BLT (Bead -Linked Transposome)** | 10 µL per sample   | 2&deg;C to 8&deg;C     |
| **TSB (Tagment Stop Buffer)**      | 10 µL per sample   | 15&deg;C to 30&deg;C   |
| **TWB (Tagment Wash Buffer)**      | 300 µL per sample  | 15&deg;C to 30&deg;C   |
| **1.7 mL microcentrifuge tubes**   | Varies             | Room temperature       |
| **96-well 0.2 mL PCR plate**       | Varies             | Room temperature       |
| **Microseal "A" and "B" film**     | Varies             | Room temperature       |

**Preparation**

- If not completed already, **carefully** pool an equal amount of each gene's corresponding sample into the same well of a half skirt plate to reach a total of 30 µL before continuing

- Add appropriate volume of DNA as obtained from gDNA dilution calculator to each well of a 96-well PCR plate so that the total input amount (100 - 1000 ng) is within the desired range. If DNA volume <30 µL, add nuclease-free molecular grade water to the DNA samples to bring the total volume to 30 µL. Add DNA to the molecular-grade water and mix well by gently pipetting approximately 5-10 times

- Vortex BLT for a minimum of 10 seconds and ensure proper suspension of beads, repeat if necessary. Do not centrifuge

- Be sure all samples are mixed thoroughly by pipetting gently 10 times.

- Bring BLT (from refrigerator) and TB1 (from freezer) to room temperature.

- **NOTE:** _Ensure that BLT is stored upright at all times, so that the beads remain submerged in the buffer. BLT must be stored at temperatures above 2°C. Do not use BLT that has been stored below 2°C and ensure that it is never frozen._

- Label a 96-well PCR plate, or equivalent, with Run ID

- Set up thermal cycler and choose preheat lid option:

**Thermocycler Program: "Flex 1" with reaction volume set to 50ul and choose the preheat lid option set to 100°C**

| Flex 1                 |
| :---:                  |
|55&deg;C for 15 min     |
|10&deg;C for Infinity   |

- Check TSB for precipitate (if present, warm at 37°C for up to 10 minutes and vortex) and ensure it is at room temperature

- Set up thermal cycler and choose preheat lid option:

**Thermocycler Program: "Flex 2" with reaction volume set to 60ul and choose the preheat lid option set to 100°C**

| Flex 2                 |
| :---:                  |
|37&deg;C for 15 min     |
|10&deg;C for Infinity   |


**Table 9a. PCR Tagmentation Master Mix Table 9a**

| **Reagent** | **Volume (µL) per sample** |
| ---         | ---                        |
| TB1         | 10 µL                      |
| BLT         | 10 µL                      |


**Procedure: Tagment Genomic DNA**

1. Prepare tagmentation master mix based on Table 9a above.

2. Vortex the tagmentation master mix well. If the master mix was set up in a trough, pipette mix instead.

3. Add 20 µL of master mix to each sample well and mix well by resuspending the beads 10 times. **Do not spin the plate.**

5. Seal the plate with Microseal B (or equivalent) and incubate the plate on the pre-programmed thermal cycler setting "Flex 1" with volume set to 50 µL and lid heated option at 100°C:

| Flex 1                 |
| :---:                  |
|55&deg;C for 15 min     |
|10&deg;C for Infinity   |

- **NOTE:** **PLEASE PROCEED TO NEXT PROCEDURE (Post Tagmentation Clean-Up)** This is **NOT** a safe stopping point in the procedure and post tagmentation clean up should be commenced once the samples have reached 10°C.

<a id="tag_cleanup"></a>

#### II: Post Tagmentation Clean-Up ####

1. Again, check TSB for precipitate (if present, warm at 37°C for up to 10 minutes and vortex) and ensure it is at room temperature.

2. Add 10 µL of TSB to each sample with a multichannel pipet and pipet gently 10 times to mix and re-suspend the beads.

4. Seal the plate with Microseal A (or equivalent) and incubate the plate on the pre-programmed thermal cycler setting "Flex 2" with volume set to 60 µL and lid heated option at 100°C:

| Flex 2                 |
| :---:                  |
|37&deg;C for 15 min     |
|10&deg;C for Infinity   |

- **NOTE:** **PLEASE PROCEED TO NEXT STEPS** This is **not** a safe stopping point in the procedure and it is recommended to proceed to the washing steps after samples have reached 10°C.

4. While samples are incubating, thaw EPM on ice and thaw indices at room temperature.

5. After incubation, remove from thermal cycler, quick spin the plate, remove microseal, and transfer the 60 µL sample volumes to a new 96‐well U-Shaped-Bottom Microplate, and place on a magnet for 3 minutes until solution is clear (or until beads form a tight pellet).

6. Using a multichannel pipette set at 100 µL, remove and discard supernatant.

7. Complete steps **A** - **D** below two times:

    - A. Remove the sample plate from the magnetic stand and add 100 µL TWB directly onto the beads.

    - B. Set multichannel pipet to 90 µL and pipette slowly until beads are fully resuspended. If necessary, scrape the side of the well with the pipette tips to re-suspend the beads.

    - C. Place the plate on the magnetic stand and wait until the solution is clear (~3 minutes).

    - D. Using a multichannel pipette set at 110 µL, remove and discard supernatant.

- **NOTE:** _To minimize the potential of TWB foaming during the tagmentation wash, use a deliberately slow pipetting technique to avoid incorrect volume aspirauion and incomplete mixing._

8. Remove the plate from the magnetic stand and add 100 µL TWB.

9. Pipette each sample well slowly to resuspend the beads.

10. Place on the magnetic stand until the solution is clear (~3 minutes). Allow TWB to remain in the wells (to prevent drying of beads) and proceed to amplification steps.

**PLEASE PROCEED TO NEXT PROCEDURE (Amplification of Tagmented DNA (Library Indexing)). <ins>THIS IS NOT A RECOMMENDED SAFE STOPPING POINT.</ins>**


<a id="lib_indexing"></a>

#### III: Amplification of Tagmented DNA (Library Indexing) ####

This step amplifies the tagmented DNA using a limited-cycle PCR program. The PCR step adds Index 1 (i7) adapters, Index 2 (i5) adapters, and sequences required for sequencing cluster generation.

**Consumables**

**Table 10. Amplification of Tagmented DNA (Library Indexing) Consumables**

| Item                            | Quantity         | Storage               |
| ---                             | ---              | ---                   |
| **EPM (Enhanced PCR Mix)**      | 20 µL per sample | -15&deg;C to -25&deg;C|
| **Index 1 adapters (plate)**    | 10 µL per sample | -15&deg;C to -25&deg;C|
| **Nuclease-free water**         | 20 µL per sample | Room temperature      |
| **1.7 mL microcentrifuge tube** | Varies           | Room temperature      |
| **Microseal "A" film**          | Varies           | Room temperature      |

> Index 1 adaptors: Catalog #20018708 (96 samples) or #20018707 (24 samples)

**Preparation**

- Using the samples suspended in the 100 µL of TWB from the Post Tagmentation Clean-Up step, follow the procedure outlined below.

- Set up thermal cycler and choose preheat lid option:

**Thermocycler Program: "Flex 3" with reaction volume set to 50ul and choose the preheat lid option set to 100°C:**

| Step   | Temperature                                 | Time (min) |
| :--:   |:----:                                       | :----:     |
| 1      | 68&deg;C                                    | 3:00       |  
| 2      | 98&deg;C                                    | 3:00       |  
|        |                                             |            |
| **3**  | **98&deg;C**                                | **0:45**   |      
| **4**  | **62&deg;C**                                | **0:30**   |      
| **5**  | **68&deg;C**                                | **2:00**   |    
|        | **Repeat Steps 2-4 for 4 cycles</br> (5 total)** |            |  
|        |                                             |            |
| 6      | 68&deg;C                                    | 1:00       |
| 7      | 4&deg;C                                     | Infinity   |

**Table 10. PCR Master Mix**

| **Reagent**           | **Volume (µL) per sample** |
| ---                   | ---                        |
| EPM                   | 20 µL                      |
| Molecular Grade Water | 20 µL                      |

**Procedure**

1. Briefly vortex the thawed EPM immediately before use.

2. Prepare the PCR master mix based on **Table 10** above.

    **NOTE:** _It is recommended to increase the number of samples during master mix calculation by 10% to ensure sufficient master mix volume_.

3. Vortex and quick spin the PCR master mix.

4. Using a multichannel pipette set at 200 µL remove TWB from beads. Use a small volume pipette to ensure removal of residual TWB before proceeding.

    **NOTE:** _Removal of TWB is crucial, as it can impede PCR. However, any foam remaining on the wells will not negatively impact the library._

5. Remove from the magnet and immediately add 40 µL of PCR master mix to each sample and gently pipet to mix, re-suspending the pellet. If necessary, scrape the side of the well with the pipette tips to resuspend the beads. Transfer the 40 µL sample volumes to a new PCR plate.

6. Add 10 µL of appropriate index pair from indices plate to each sample well.

    **NOTE:** It is recommended to pierce the foil of the desired well on the index plate with a new 200 µL pipet tip, then to use a fresh pipette tip to withdraw the indices from the wells, followed by re- sealing the index plate with a new foil cover (i.e. Microseal F) after each use. Make sure that the index is oriented correctly. Handle plate gently to maintain index at the bottom of the plate. If not, spin plate to make sure that index is towards bottom of the plate.

     **NOTE:** Index should be added as next available down the columns.

7. Using a multichannel pipette set at 40 µL mix by pipetting a minimum of 10 times.

8. Seal the plate with Microseal A (or equivalent) and place the plate on the pre-programmed thermal cycler setting "Flex 3" with volume set to 50 µL and lid heated option at 100°C.

**SAFE STOPPING POINT**

The plate may be sealed with Microseal B or equivalent and stored at 2°C to 8°C for up to 3 days. If you choose to continue, please proceed to Library PCR Clean-Up.

 * * *

<a id="lib_cleanup"></a>

#### IV: NGS Library Clean‐Up ####

This step uses Sample purification beads to clean up the final library before quantification.

**Consumables**

**Table 11. Library PCR Purification Consumables**

| Item | Quantity | Storage |
| --- | --- | --- |
| **RSB (Resuspension Buffer)** | 32 µL per sample | -15° to -25°C (after initial thaw, can keep at 2° to 8°C |
| **SPB (Sample Purification Beads)** | 40.8 µL per sample | 2° to 8°C |
| **Freshly Prepared 80% Ethanol (EtOH)** | 44.2 µL per sample | Room temperature |
| **96‐well 0.2 mL PCR plate** | Varies | Room temperature |
| **Nuclease-free water** | Varies | Room temperature |
| **Microseal 'B' film and 'F' foil** | Varies | Room temperature |
| **96‐well U-Shaped-Bottom Microplate** | Varies | Room temperature |
| **96-well 0.2 mL PCR plate** | Varies | Room temperature |
| **1.7 mL microcentrifuge tube** | Varies | Room temperature |

**Preparation**

- Thaw RSB at room temperature and vortex to mix. For future use, it is recommended to aliquot RSB into a smaller volume container (such as a 1.5 mL tube) to reduce thawing times.

- Bring the SPB to room temperature- wait at least 30 minutes. Once at room temperature, vortex and invert the beads several times to re-suspend any particles that may have settled. Ensure magnetic beads are well (evenly) distributed before adding them to samples.

- **NOTE:** Proceed very cautiously during this procedure and take your time to ensure as little bubble formation as possible.

-  Prepare a fresh dilution stock of 80% ethanol sufficient for all samples:

| **Reagent**                | **Volume (ml) per sample** | **Example: 20 samples** |
| ---                        | ---                        | ---                     |
| **100% ethanol**           | 0.4                        | 8 mL                     |
| **Molecular grade water**  | 0.1                        | 2 mL                     |

- **NOTE:** 80% ethanol is hygroscopic. That is, when opened the ethanol will both evaporate and absorb water over time. Re-use eventually will be at a lower concentration. There is also miscibility involved with ethanol and water. For example, measuring out 80 mL of ethanol and topping off to 100 mL with water will generate ~65% ethanol. Measuring 80 mL ethanol and 20 mL water separately, then combining them will generate ~95 mL of 80% ethanol. Make sure to use molecular biology grade water (DNAse, RNase and Protease free).

**Table 11. SPB Master Mix**

| **Reagent**               | **Volume (µL) per sample** |
| ---                       | ---                        |
| **SPB**                   | 40.8                       |
| **Molecular grade water** | 44.2                       |

- **NOTE:** _It is recommended to increase the number of samples by 10% to ensure sufficient volume of master mix._


**Procedure**

1. Centrifuge the Library Amplification and Index PCR plate at 280 × g at 20°C for 1 minute to collect condensation, carefully remove seal.

2. Prepare SPB master mix in a 2 mL tube based on **Table 11a** above.

3. Transfer the 50 µL sample volumes to a new 96‐well U-Shaped-Bottom Microplate. Place sample plate on the magnet for 5 minutes (or until beads have formed a tight pellet).

4. Transfer 45 µL of supernatant (now containing the DNA) to new deep well plate.

5. Remove sample plate from the magnet.

6. Vortex SPB master mix thoroughly.

7. Using a multichannel pipette mix briefly and add 85 µL of SPB master mix to each sample and pipette mix a minimum of 10 times.

   **NOTE:** _Use caution when mixing as the volume will be >100 µL._

8. Incubate at room temperature for 5 minutes **.**

9. Place on the magnet for 3-5 minutes (or until beads form a tight pellet).

10. During incubation re-vortex the stock SPB.

11. After incubation, while keeping the plate still on the magnet, transfer 120 µL of supernatant - which now contains the library DNA - to new wells. _If desired, up to 125 µL of supernatant can be transferred_

12. Remove the plate from the magnet and add 14.4 µL of stock SPB solution to the supernatant. If volume other than 105 µL was used, then maintain a bead ratio of 0.12x.

13. With multichannel pipette set to 100 µL, gently pipet 10 times to mix.

14. Incubate at room temperature for 5 minutes.

15. Place on magnet for 3-5 minutes (or until beads form a tight pellet and supernatant clears).

16. With multichannel pipette set to 200 µL, remove and discard supernatant (DNA is now bound to the beads).

17. With the Library amplification plate on the magnetic stand, perform steps **A** - **C** below **twice (for a total of two washes)**:

    A. Add 170 µL of fresh 80% ethanol. (DO NOT add directly to the bead, and DO NOT mix).

    B. Incubate the plate on the magnetic stand for 30 seconds.

    C. Carefully remove and discard all the ethanol.

18. Use a P20 multichannel pipette with fine pipette tips to remove excess ethanol.

19. With the Library amplification plate still on the magnetic stand, allow the beads to air‐dry for 3-5 minutes.

    **NOTE**: Make sure not to over dry the beads. Bead pellets will appear cracked if over dried. If cracking is observed, immediately re-suspend beads as described below regardless of drying time.

20. Remove the plate from the magnetic stand and add 32 µL RSB to each well of the plate and gently pipet a minimum of 10 times to thoroughly mix.

21. Incubate at room temperature for 2-5 minutes.

22. Place the plate back on the magnetic stand for 3 minutes or until the supernatant has cleared.

23. Using a multichannel pipette, carefully transfer 25 µL of the supernatant from the Library amplification plate to a new 96‐well PCR plate. Change tips between samples to avoid cross‐contamination.

**SAFE STOPPING POINT**
> If you do plan to stop here, seal the plate with Microseal "B" adhesive seal. Store the plate at ‐15° to ‐25°C for up to a week.


 * * *

<a id="lib_clustering"></a>

## Library Clustering ##

#### [Part I: Library Pooling](#lib_pool) ####
#### [Part II: Library Quantification](#lib_quant) ####
#### [Part III: Library Normalization](#lib_norm) ####

It is important to consider library size when preparing samples for cluster generation. Because the clustering process preferentially amplifies shorter libraries in a mixture of fragments, large libraries tend to cluster less efficiently than small libraries. The DNA concentration used for clustering can be adjusted to increase the cluster density of larger libraries. Consider table 1 below:

**Library Denaturing and MiSeq Sample Loading**

|Average Library Size | Conversion Factor    |  DNA Concentration for Cluster Generation  |
|---                  | --                   |--                                          |
| 250 bp              |  1 ng/µL = 6.0 nM    |      6 - 12 pM                             |
| 500 bp              |  1 ng/µL = 3.0 nM    |      6 - 12 pM                             |
| 1,000 bp - 1,500 bp |  1 ng/µL = 1.5 nM    |      6 - 12 pM                             |

>The values presented here are approximations, and exact vallues determined for each experiiment may differ from these guidelines.


**Procedure**

<a id="lib_pool"></a>

#### Part I: Library Pooling ####

1. Aliquot 5 µL of diluted DNA from each library into a 1.5 microcentrifuge tube and mix aliquots for pooling libraries with unique indices. Depending on coverage needs, up to 384 libraries can be pooled for one MiSeq run.


<a id="lib_quant"></a>

#### Part II: Library Quantification ####

**Background**
Quantification of fragment size and concentration to determine library concentration in nM. Illumina recommends quantifying your libraries using a fluorometric quantification method that uses dsDNA binding dyes.

- Our laboratory uses the Agilent D5000 ScreenTape System Quick Guide protocol from Agilent Technologies to determine the fragment size of our libraries.
- Our laboratory uses the Qubit® dsDNA HS Assay Kits protocol from Life Technologies to determine the conentration of our libraries.

**DNA Concentration in nM**

After determining the fragment size and concentration of your pooled product, you will calculate the DNA concentration in nM, based on the size of DNA amplicons as determined by an Agilent Technologies 2100 Bioanalyzer trace and concentration by Qubit as follows:

(concentration in ng/µL)(10^6) / (660 g/mol)(average library size) = concentration in nM

For example:
(15 ng/µL)(10^6) / (660 g/mol)(500 bp) = 45 nM

**Agilent Technologies Agilent D5000 ScreenTape System**

This SOP format was adapted from the Agilent D5000 ScreenTape System Quick Guide [protocol](https://www.agilent.com/cs/library/usermanuals/public/HS-D5000_QuickGuide.pdf) from Agilent Technologies.

**Consumables**

**Table 12. TapeStation Consumables**

| Item              | Quantity                  | Storage   |
| ---               | ---                       | ---       |
| **Sample Buffer** | 10 µL per sample          | 2° to 8°C |
| **D5000 Ladder**  | 1 µL                      | 2° to 8°C |
| **ScreenTape**    | Holds 16 samples per tape | 2° to 8°C |

**Prepare TapeStation System D5000**

- Launch the 4150 TapeStation Controller Software.

- Load single D5000 ScreenTape device and loading tips into the 4150 TapeStation instrument.

**Procedure:** _Determine fragment size_

**Sample Preparation D5000 ScreenTape Assay**

1. Allow reagents to reach equilibrium at room temperature for 30 minutes.

2. Vortex mix before use.

3. Prepare ladder by mixing 2 µL D5000 Sample Buffer (green lid) with 2 µL D5000 Ladder (yellow lid) in a tube strip.

4. Prepare sample by mixing 10 µL D5000 Sample Buffer (green lid) with 1 µL DNA sample in different tube strips.

5. Spin down, then vortex using IKA vortexer and adaptor at 2000 rpm for 1 minute.

6. Spin down to position the sample at the bottom of the tube.

**Sample Fragment Size Analysis**

1. Load samples into the 4150 TapeStation instrument.

2. Select the required samples on the 4150 TapeStation Controller Software.

3. Click Start and specify a filename with which to save your results.

**SAFE STOPPING POINT**
> If you do not plan to proceed to _Part II Qubit Flurometer 3.0 dsDNA HS Assay_, leave your sample in 4°C for maximum of a week.


**Qubit Fluorometer 3.0 dsDNA HS Assay**

This SOP format was adapted from the Qubit® dsDNA HS Assay Kits [protocol](https://tools.thermofisher.com/content/sfs/manuals/Qubit_dsDNA_HS_Assay_UG.pdf) from Life Technologies.

**Consumables:**

**Table 13. Qubit 3.0 Fluorometer Consumables**

| Item                      | Quantity                               | Storage          |
| ---                       | ---                                    | ---              |
| **Qubit dsDNA HS Buffer** | 199 µL per sample for working solution | Room temperature |
| **Qubit dsDNA HS Reagent**| 1 µL per 199 µL of HS Buffer           | Room temperature |
| **Standard #1**           | 10 µL per use                          | 2° to 8°C        |
| **Standard #2**           | 10 µL per use                          | 2° to 8°C        |
| **Qubit™ Assay Tubes**    | 1 per sample and 1 for each ladder     | Room temperature |

**Before you begin**

- The final volume in each tube must be 200 µL.

- Each standard tube requires 190 µL of Qubit working solution + 10 µL of the standard.

- Each sample tube requires anywhere from 180–199 µL + the corresponding volume to complete the necessary 200 µL.
  - This laboratory uses 195 µL working solution + 5 µL of sample.

- Careful pipetting is critical to ensure that the exact volume is added to the working solution—work SLOWLY.

- Be sure to use a clean plastic tube each time you prepare Qubit working solution. Do not mix the working solution in a glass container.

**Procedure:** _Determine library concentration_

Standard and Sample Preparation

1. Prepare the tubes: Set up two 0.5-mL tubes for standards, and the required number of tubes for samples.

- **Note** Use only the thin-wall, clear, 0.5-mL PCR tubes (described in **Table 2** User‐Supplied Consumables).

2. Label the tube lids. _Do not label the side of the tube as this could interfere with the sample read._

3. Prepare sufficient Qubit working solution to accommodate all standards and samples by diluting the Qubit dsDNA HS Reagent 1:200 in Qubit dsDNA HS Buffer. 1 µL Qubit dsDNA HS Reagent + 199 µL Qubit dsDNA HS Buffer. For example, for 8 samples, prepare enough working solution for the samples and two standards: ~200 µL per tube in 10 tubes yields 2 mL of working solution (10 µL of Qubit reagent plus 1990 µL of Qubit buffer).

3. Prepare the standards by adding 190 µL of Qubit working solution to each of the tubes used for standards. Add 10 µL of each Qubit standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.

4. Prepare the samples by adding Qubit working solution to individual assay tubes so that the final volume in each tube after adding the sample is 200 µL.

- **NOTE:** Your sample added should be between 1–20 µL and Qubit working solution between 180–199 µL (**See Table 13a** below).

5. Add each sample to the assay tubes containing the correct volume of Qubit working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 µL.

6. Allow all tubes to incubate at room temperature for 2 minutes.

**Table 13: Sample and Quibit Working Solution**  

| Working Solution Volume | 199 µL | 195 µL | 190 µL | 185 µL | 180 µL |
| ---                     | ---    | ---    | ---    | ---    | ---    |
| Sample Volume           | 1 µL   | 5 µL   | 10 µL  | 15 µL  | 20 µL  |


**Standard and Sample Reading**

1. On the home screen of the Qubit 3.0 Fluorometer, select "**dsDNA**", then "**High Sensitivity**", and then "**Read Standards**."

- **NOTE:** If you have already performed a calibration for the selected assay, the instrument prompts you to choose between reading new standards and running samples using the previous calibration. If you wish to use the previous calibration, disregard step 3 in the Standard and Sample Preparation step, and skip to step 4 below. Otherwise, continue with step 2 below.

2. Insert the tube containing Standard #1 into the sample chamber, close the lid, and then press **Read standard**. When the reading is complete (~3 seconds), remove Standard #1.

3. Insert the tube containing Standard #2 into the sample chamber, close the lid, and then press **Read standard**. When the reading is complete, remove Standard #2.

4. Press **Run samples**.

5. On the assay screen, select the sample volume and units using the + or – buttons on the wheel to select the sample volume added to the assay tube (from 1–20 µL).

6. From the dropdown menu, select the units for the output sample concentration (ng/µL).

7. Insert a sample tube into the sample chamber, close the lid, and press Read tube. When the reading is complete (~3 seconds), remove the sample tube and repeat until all samples have been read.


<a id="lib_norm"></a>

#### Part III: Library Normalization ####

Dilute concentrated final library using Resuspension Buffer (RSB) or fresh 10 mM Tris pH 8.5 to **4 nM.**

For example:
Given a calculated concentration of 45nM, use (C1)(V1) = (C2)(V2) to calculate how much RSB and sample to mix to create a 4nM concentration:

45nM(V1) = 4nM(20 µL)
V1 = 1.78 µL of sample + 18.22 µL of RSB = 20ul of a 4nM concentration

**SAFE STOPPING POINT**
>If you do not plan to proceed to _Library Denaturing and MiSeq Sample Loading_, leave your sample in 4°C for a maximum of one week.


 * * *

## Sample sequencing ##

<a id="sample_loading"></a>

#### Library Denaturing and MiSeq Sample Loading ####

In preparation for cluster generation and sequencing, pooled libraries are denatured with NaOH, diluted with hybridization buffer, and then heat denatured before MiSeq sequencing. Each run must include a minimum of 5% PhiX to serve as an internal control for these low-diversity libraries. Illumina recommends using MiSeq v2 reagent kits for improved run metrics.

**Consumables**

**Table 14. Library Denaturing and Miseq Sample Loading Consumables**

| Item                                  | Quantity         | Storage          |
| ---                                   | ---              | ---              |
| **RSB (Resuspension Buffer)**         | 6 µL             | -15° to -25°C    |
| **HT1 (Hybridization Buffer)**        | 1540 µL          | -15° to -25°C    |
| **0.2 N NaOH (less than a week old)** | 10 µL            | Room temperature |
| **200 mM Tris-HCl pH7.0**              | 5 µL             | Room temperature |
| **PhiX Control Kit v3 (FC‐110‐3001)** | 2 µL             | -15° to -25°C    |
| **MiSeq v2 reagent cartridge**        | 1 cartridge      | -15° to -25°C    |
| **1.7 mL microcentrifuge tubes (screw cap recommended)** | 3 tubes          |
| **2.5 L ice bucket**                  |                  |                  |

**Preparation**

1. Begin thawing the reagent cartridge and HT1 before denaturing and diluting libraries by placing them in a room temperature water bath for about an hour.

2. Once thawed, store the cartridge and HT1 in the ice bucket until ready for sample loading.

3. Obtain an ice bucket for your thawed cartridge, freshly made reagents, and sample.

4. Check pH of the stock 1.0N NaOH and the resulting 0.2N NaOH dilution using pH reader.

- **NOTE:** CO2 in the room will acidify NaOH over time. It is absolutely **critical** that the NaOH has a minimum pH of 12.5.

5. Prepare a fresh dilution of 0.2 N NaOH.

- **NOTE: this is a critical step; NaOH must be prepared fresh every time.**

6. Using a 1000 µL pipette, measure out 800 µL of laboratory-grade water.

7. In a separate microcentrifuge tube, measure 200 µL of stock 1.0N NaOH.

8. Combine the two volumes and then invert several times to mix.

- **NOTE:** This results in a 1 mL of 0.2N NaOH solution; the resulting solution cannot be stored and must be used within 6 hours.

- **NOTE:** The final concentration of NaOH cannot exceed 0.001 (1 mM). _Higher NaOH concentrations will inhibit library hybridization to the flow cell and result in very low cluster density._

9. If you have not already done so, prepare a 200 mM stock of Tris-HCl pH7.0 by combining 800 µL of Laboratory-grade water and 200 µL of Tris-HCl 1M.


**Denature DNA**

1. Combine the 4nM pooled library (5 µL) and 0.2N NaOH (5 µL) in a microcentrifuge tube:

2. Set aside the remaining dilution of 0.2 N NaOH to prepare a PhiX control within the next 12 hours.

3. Vortex briefly to mix the sample solution, and then centrifuge the sample solution at 280 × g (or about 1500rpm) at 20°C for 1 minute.

4. Incubate for 5 minutes at room temperature to denature the DNA into single strands.

5. To the 10 µL of denatured library, add 5 µL of 200 mM Tris-HCl pH7.0 to neutralize the NaOH.

6.  Add pre-chilled HT1 (985 µL) to the denatured DNA + Tris-HCl (15 µL). Adding the HT1 results in a 20 pM denatured library in 1 mM NaOH.

7. Place the denatured DNA on ice until you are ready to proceed to final dilution.


|    Quick Review/Guide for denaturing 4nM library                             |                                    
| :---:                                                                        |
| **_Step 1_**                                                                 |  
| **5 µL of 4 nM library + 5 µL of 0.2N NaOH + 5 µL of 200 mM Tris-HCl pH 7.0** |
|  **=**                                                                       |
| 15 µL of 1.30 nM library 0.067N NaOH + 66.7 mM Tris-HCl pH 7.0                |
| **_Step 2_**                                                                 |  
|  **Add 985 µL of chilled HT1**                                               |  
|  **=**                                                                       |
| 1 mL of **0.001N NaOH** and 20 pM denatured library + 1 mM Tris-HCl pH 7.0      |


**NOTE: If you have to start with a lower concentration library, follow the below protocol for denaturing a 2nM library.**


|    Quick Review/Guide for denaturing 2nM library                             |                                    
| :---:                                                                        |
| **_Step 1_**                                                                 |  
| **5 µL of 2 nM library + 5 µL of 0.2N NaOH + 5 µL of 200 mM Tris-HCl pH 7.0** |
|  **=**                                                                       |
| 15 µL of 0.67 nM library 0.067N NaOH + 66.7 mM Tris-HCl pH 7.0                |
|  **_Step 2_**                                                                |
|  **Add 985 µL of chilled HT1**                                               |  
|  **=**                                                                       |
| 1 mL of **0.0005N NaOH** and 10pM denatured library + 1 mM Tris-HCl pH 7.0     |


**Dilution chart for 10pM library:**

| **Final Concentration**     | 6pM     | 8pM    |10pM    |
|---                          |---      |---     |---     |
| **10pM denatured library**  | 360 µL  |	480 µL | 600 µL	|
| **Pre-chilled HTI** |240 µL |	120 µL  |	0 µL   |        |

**Dilute Denatured DNA**

 - **NOTE:** Illumina recommends targeting 1000–1200 K/mm² raw cluster densities using MiSeq v2 reagents.

1. Dilute the denatured DNA to the desired concentration using the following example:

2. Invert several times to mix and then pulse centrifuge the DNA solution.

3. Place the denatured and diluted DNA on ice.

**Denature and Dilute of PhiX Control**

Use the following instructions to denature and dilute the 10 nM PhiX library to the same loading concentration as the Amplicon library. The final library mixture must contain at least 5% PhiX.

1. Combine 10 nM PhiX library (2 µL) and RSB (3 µL) to dilute the PhiX library to 4 nM:

2. Combine 4 nM PhiX library (5 µL) and 0.2 N NaOH (5 µL) of 4 nM PhiX and 0.2 N NaOH in a microcentrifuge tube.

3. Vortex briefly to mix the 2 nM PhiX library solution.

4. Incubate for 5 minutes at room temperature to denature the PhiX library into single strands.

5. To the 10 µL of denatured library, add 5 µL of 200 mM Tris-HCl pH7.0 to neutralize the NaOH.

6. Add denatured PhiX library (15 µL) and Pre‐chilled HT1 (985 µL) to make a 20 pM PhiX library:

7. Dilute the denatured 20 pM PhiX library to the same loading concentration as the Amplicon library as shown in **Table 14a** below.  

8. Invert several times to mix and then pulse centrifuge the DNA solution.

9. Place the denatured and diluted PhiX on ice.

**Table 14a: Clustering Library Table**

| **Final Concentration**     | 2pM     | 4pM    | 6 pM   | 8 pM   |
|---                          |---      |---     |---     |---     |
| **10pM denatured library**  | 60 µL  |	120 µL | 180 µL	| 240 µL |
| **Pre-chilled HTI**         | 540 µL |	480 µL | 420 µL | 360 µL |    

**Table 14b: Clustering Library Table continued**

| **Final Concentration**     | 10pM     | 12pM    | 15pM   | 20pM   |
|---                          |---       |---      |---     |---     |
| **10pM denatured library**  | 300 µL   |	360 µL | 450 µL	| 600 µL |
| **Pre-chilled HTI**         | 300 µL   |	240 µL | 150 µL | 0 µL   |  



**Combine Amplicon Library and PhiX Control**

- **NOTE:** The recommended PhiX control spike‐in of greater than 5% for low diversity libraries is possible with RTA v1.17.28 or later, which is bundled with MCS v2.2. For optimal performance, update to v3 software (MCS 2.3). If you are using an older version of the MiSeq software or sequencing these libraries on the GA or HiSeq, Illumina recommends using greater than 25% PhiX control spike‐in.

1. Combine denatured and diluted PhiX control (30 µL) and denatured and diluted amplicon library (570 µL). This will result in a 5% PhiX spike-in.  

2. Set the combined sample library and PhiX control aside on ice until you are ready to load the mixture into the MiSeq v2 reagent cartridge.

3. Invert the tube 1–2 times to mix and load all 600ul into the designated well in the cartridge.


<a id="supporting_info"></a>

## Supporting Information ##

The protocols described in this guide assume that you are familiar with the contents of this section and have obtained all of the requisite equipment and consumables.

**Acronyms**

**Table 15. Definitions and Acronyms**

| Acronym | Definition |
| --- | --- |
| **PCR** | Polymerase Chain Reaction- a technique used to amplify 1 to a few copies of a piece of DNA across several orders of magnitude, generating thousands to millions of copies of a single DNA strand |
| **Primer** | A strand of short nucleic acid sequences that serves as a starting point for DNA synthesis during PCR |
| **Amplicon** | A piece of amplified DNA that is the product of a PCR reaction |

## Illumina DNA/RNA UD Indexes

Each sample processed during library preparation will be identified post-sequencing by the unique index that is molecularly attached to that sample during the [**Amplification of Tagmented DNA (Library Indexing)**](#lib_indexing) steps. The indexes used in this protocol are the [IDT for Illumina DNA/RNA UD Indexes](https://support.illumina.com/sequencing/sequencing_kits/idt-nextera-dna-udi/product-files.html) and for each sequence of the 384 indexes there is a [Kit Definition File](https://support.illumina.com/documents/downloads/productfiles/unique-dual-indexes/nextera-dna-udi/nextera-dna-udi-lrm-library-kit-definition-iSeq-MiSeq-nextera-dna-flex-set-a-b-c-d-2x151-384-samples.zip) and [Sample Sheet File](https://support.illumina.com/documents/downloads/productfiles/unique-dual-indexes/nextera-dna-udi/nextera-dna-udi-lrm-samplesheet-iSeq-MiSeq-nextera-dna-flex-set-a-b-c-d-2x151-384-samples.zip) available.


## Prevent PCR Contamination ##

The PCR process is commonly used in the laboratory to amplify specific DNA sequences. Unless proper laboratory hygiene is used, PCR products can contaminate reagents, instrumentation, and genomic DNA samples, causing inaccurate and unreliable results. PCR product contamination can shut down lab processes and significantly delay normal operations.

Make sure that the lab is set up appropriately to reduce the risk of PCR product contamination:

* Physically Separate Pre-PCR and Post-PCR Areas
  * Physically separate laboratory space where pre‐PCR processes are performed (DNA extraction, quantification, and normalization) from the laboratory space where PCR products are made and processed (post‐PCR processes).
  *  Never use the same sink to wash pre‐PCR and post‐PCR troughs.<
  * Never share water purification systems for pre‐PCR and post‐PCR processes.
  * Store all supplies used in the protocols in the pre‐PCR area, and transfer to the post‐ PCR area as needed.
    Use Dedicated Equipment and Supplies
  * Dedicate separate full sets of equipment and supplies (pipettes, centrifuges, oven, heat block, etc.) to pre‐PCR and post‐PCR lab processes, and never share between processes.
  * Dedicate separate storage areas (freezers and refrigerators) to pre‐PCR and post‐PCR consumables.

Because the pre‐ and post‐amplification reagents are shipped together, it is important to unpack the reagents in the pre‐PCR lab area. After unpacking the reagents, move the post-amplification reagents to the proper post‐PCR storage area.

**Pre‐PCR and Post‐PCR Lab Procedures**

To prevent PCR product contamination, it is important to establish lab procedures and follow best practices. Illumina recommends daily and weekly cleaning of lab areas using 0.5% Sodium Hypochlorite (10% Bleach).

  **CAUTION:**  To prevent sample or reagent degradation, make sure that all vapors from the cleaning solution have fully dissipated before beginning any processes.

**Daily Cleaning of Pre‐PCR Area**

A daily cleaning of the pre‐PCR area using a 0.5% Sodium Hypochlorite (10% Bleach) solution helps to eliminate PCR product that has entered the pre‐PCR area. Identify pre‐PCR areas that pose the highest risk of contamination, and clean these areas with a 0.5% Sodium Hypochlorite (10% Bleach) solution before beginning any pre‐PCR processes. High‐risk areas might include, but are not limited to, the following items:
  * Benchtops
  * Door handles
  * Refrigerator/freezer door handles
  * Computer mouse
  * Keyboards

**Daily Cleaning of Post‐PCR Area**

Reducing the amount of PCR product in the post‐PCR area helps reduce the risk of contamination in the pre‐PCR area. Daily cleaning of the post‐PCR area using a 0.5% Sodium Hypochlorite (10% Bleach) solution helps reduce the risk of contamination. Identify post‐PCR areas that pose the highest risk of contamination, and clean these areas with a 0.5% Sodium Hypochlorite (10% Bleach) solution daily. High‐risk areas might include, but are not limited to, the following items:

  * Thermal cyclers
  * Bench space used to process amplified DNA
  * Door handles
  * Refrigerator/freezer door handles
  * Computer mouse
  * Keyboards

**Weekly Cleaning of All Lab Areas**

One time a week, perform a thorough cleaning of the pre‐PCR and post‐PCR areas using 0.5% Sodium Hypochlorite (10% Bleach).

  * Clean all benchtops and laboratory surfaces.
  * Clean all instruments that are not cleaned daily.
  * Thoroughly mop lab floors.
  * Make sure that personnel responsible for weekly cleaning are properly trained on prevention of PCR product contamination.

**Items Fallen to the Floor**

The floor is contaminated with PCR product transferred on the shoes of individuals coming from the post‐PCR area; therefore, anything falling to the floor must be treated as contaminated.

  * Disposable items that have fallen to the floor, such as empty tubes, pipette tips, gloves, lab coat hangers, must be discarded.
  * Non‐disposable items that have fallen to the floor, such as a pipette or an important sample container, must be immediately and thoroughly cleaned. Use a 0.5% Sodium Hypochlorite (10% Bleach) solution to remove PCR product contamination.
  * Clean any lab surface that has come in contact with the contaminated item. Individuals handling anything that has fallen to the floor, disposable or non‐disposable, must discard their lab gloves and put on a new pair.

## Best Practices ##

When preparing libraries for sequencing, always adhere to good molecular biology practices. Read through the entire protocol before starting to make sure that all of the required materials are available and your equipment is programmed and ready to use.

**Handling Liquids**<br>
Good liquid handling measures are essential, particularly when quantifying libraries or diluting concentrated libraries for making clusters.

  * Small differences in volumes (±0.5 µL) can sometimes cause large differences in cluster numbers (~100,000).
  * Small volume pipetting can be a source of potential error in protocols requiring the generation of standard curves, such as qPCR, or small but precise volumes, such as the Agilent Bioanalyzer.
  * If small volumes are unavoidable, use due diligence to make sure that pipettes are correctly calibrated.
  * Make sure that pipettes are not used at the volume extremes of their performance specifications.
  * Prepare the reagents for multiple samples simultaneously, to minimize pipetting errors, especially with small volume enzyme additions. As a result, pipette one time from the reagent tubes with a larger volume, rather than many times with small volumes. Aliquot to individual samples in a single pipetting movement to allow for standardization across multiple samples.


**Handling Magnetic Beads**<br>
**NOTE:** Cleanup procedures have only been validated using the 96‐well plates and the magnetic stand specified in Tables 1 and 2. Comparable performance is not guaranteed when using a microcentrifuge tube or other formats, or other magnets.

  * Before use, allow the beads to come to room temperature.
  * Do not reuse beads. Always add fresh beads when performing these procedures.
  * Immediately before use, vortex the beads until they are well dispersed and the color of the liquid is homogeneous.
  * When pipetting beads, pipette slowly and dispense slowly due to the viscosity of the solution.
  * Take care to minimize bead loss, which can affect final yields.
  * Change the tips for each sample, unless specified otherwise.
  * Let the mixed samples incubate at room temperature for the time indicated in the protocol for maximum recovery.
  * When removing and discarding supernatant from the wells, use a single channel or multichannel pipette and take care not to disturb the beads
  * When aspirating the cleared solution from the reaction plate and wash step, it is important to keep the plate on the magnetic stand and not disturb the separated magnetic beads. Aspirate slowly to prevent the beads from sliding down the sides of the wells and into the pipette tips.
  * To prevent the carryover of beads after elution, approximately 2.5 µL of supernatant is left when the eluates are removed from the bead pellet.
  * Be sure to remove all of the ethanol from the bottom of the wells, as it can contain residual contaminants.
  * Keep the reaction plate on the magnetic stand and let it air‐dry at room temperature to prevent potential bead loss due to electrostatic forces. Allow for the complete evaporation of residual ethanol, because the presence of ethanol affects the performance of the subsequent reactions. Illumina recommends at least minutes drying time, but a longer drying time can be required. Remaining ethanol can be removed with a 10 µL pipette.
  * Avoid over drying the beads, which can impact final yields.
  * Do not scrape the beads from the edge of the well using the pipette tip.
  * To maximize sample recovery during elution, incubate the sample/bead mix for 2 minutes at room temperature before placing the samples onto the magnet.

**Avoiding Cross‐Contamination**
Practice the following to avoid cross‐contamination:

  * Open only one adapter tube at a time.
  * Change the tips for each sample, unless specified otherwise.
  * Pipette carefully to avoid spillage.
  * Clean pipettes and change gloves between handling different adapter stocks.
  * Clean work surfaces thoroughly before and after the procedure.

**Potential DNA Contaminants**

When handling and processing samples using this protocol, use best practices to avoid PCR contamination, as you would when preparing PCR amplicons.

**Temperature Considerations**

Temperature is an important consideration for making libraries:

  * Keep libraries at temperatures less than 37°C, except where specifically noted.
  * Place reagents on ice after thawing at room temperature.

**Equipment**

  * Review the programming instructions for your thermal cycler user guide to make sure that it is programmed appropriately using the heated lid function.
  * It is acceptable to use the thermal cycler tracked heating lid function.
