November 21st, 2022; version 5.3.4; updated by [Eldin Talundzic](mailto:etalundzic@cdc.gov)

> _Version 5.0: Major update made in library prep from NexteraXT to Flex protocol._

  * * *

**Table of contents**  

* [Introduction/General Overview](#intro)
* [Materials and Equipment](#materials)
* [Protocol Workflow](#workflow)
* [PET PCR Sample QC](#pet_pcr)
* [Gene PCR Enrichment](#gene_enrichment)
* [Electrophoresis](#electrophoresis)
* [PCR amplicon Clean-Up](#pcr_cleanup)
* [Tagmentation of PCR Amplicon and Tagmentation Clean Up​](#tagmentation)
* [Amplification of Tagmented DNA (Library Indexing)](#lib_indexing)
* [Library PCR Clean up](#lib_cleanup)
* [Library Pooling, Quantification, and Normalization](#lib_clustering)
* [Library Denaturing and MiSeq Sample Loading](#sample_loading)
* [Supporting Information](#supporting_info)


**IMPORTANT NOTICE**

This document provides information for an application for Illumina technology that has been demonstrated internally and may be of interest to external groups. This information is provided as‐is and is not an Illumina or CDC endorsed product and is not accompanied by any rights or warranties.

  * * *

<a id="intro"></a>

# Introduction #

**Standard Operating Procedure (SOP) describing how to prepare and sequence the full length _P. falciparum_ genes associated with antimalarial resistance on the Illumina MiSeq.**

* kelch 13 (_k13_)
* chloroquine resistant transporter (_crt_)
* multidrug resistant protein 1 (_mdr1_)
* bifunctional dihydrofolate reductase-thymidylate synthase (_dhfr_)
* dihydropteroate synthase (_dhps_)
* cytochrome b (_cytb_)

Human malaria is caused by six _Plasmodium_ species: *Plasmodium falciparum (Pf), P. vivax (Pv), P. malariae (Pm), P. ovale (Po) (P.o. curtisi and P.o. wallikeri)* and *P. knowlesi (Pk)*, which although zoonotic, can cause human infections in several South East Asian countries. Two of these, *Pf* and *Pv*, pose the greatest threat to global public health. About 3.2 billion people, half of the world's population, are at risk for malaria. In 2020, there were an estimated 241 million malaria cases in 85 malaria endemic countries, causing an estimated 627,000 deaths. In the U.S., an estimated 1,500 - 2,000 cases of malaria are imported annually. One of the greatest public health challenges for malaria control and elimination is the threat of drug resistant *Pf* parasites.

Previously effective anti-malarial treatments, chloroquine (CQ) and sulfadoxine/pyrimethamine (SP), are ineffective in many regions. Even more alarming, resistance to the least effective class of anti-malarial drugs, called artemisinins, has now emerged and spread in Southeast Asia, threatening malaria control and prevention programs globally.

Identifying and tracking drug resistance is critical for providing appropriate malaria prophylaxis and treatment guidelines. Molecular markers of resistance are available for several anti-malarial drugs, including artemisinins. Surveillance using molecular markers provides a robust system for the detection and tracking of resistant malaria parasites.

Below is a table of the major antimalarials and their associated with resistance molecular marker.

**Table 1. Antimalarials and associated resistance molecular markers**

| Antimalarial                        | Molecular Marker (Gene) | Location     |
|----                                 |----                     |----          |
| Chloroquine                         | _crt_                   | Chr 7        |
| Artemisinin                         | _k13_                   | Chr 13       |
| Atovaquone                          | _cytb_                  | Mitochondria |
| Amodiaquine, lumefantrine, quinine  | _mdr1_                  | Chr 5        |
| Pyrimethamine                       | _dhfr_                  | Chr 4        |
| Sulfadoxine                         | _dhps_                  | Chr 8        |

> Chr = chromosome; _crt+_ = chloroquine resistance transporter; _k13_ = kelch 13; _cytb_ = cytochrome b; _mdr1_ = multidrug resistance protein 1; _dhfr_ =bifunctional dihydrofolate reductase thymidylate synthase; _dhps_ = hydroxymethyldihydropterin pyrophosphokinase-dihydropteroate synthase;


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
| **1.7 ml microcentrifuge tubes** | General lab supplier |
| **10 uL barrier pipette tips** | General lab supplier |
| **10 uL multichannel pipettes** | General lab supplier |
| **10 uL single channel pipettes** | General lab supplier |
| **20 uL barrier pipette tips** | General lab supplier |
| **20 uL multichannel pipettes** | General lab supplier |
| **20 uL single channel pipettes** | General lab supplier |
| **200 uL barrier pipette tips** | General lab supplier |
| **200 uL multichannel pipettes** | General lab supplier |
| **200 uL single channel pipettes** | General lab supplier |
| **1000 uL barrier pipette tips** | General lab supplier |
| **1000 uL multichannel pipettes** | General lab supplier |
| **1000 uL single channel pipettes** | General lab supplier |
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
| **Quick-Load 1kb DNA ladder** | NEB Catalog #N0468 |
| **UltraPure™ 10X TBE Buffer** | Fisher Scientific Catalog #15581-044 |
| **AMPure XP beads for PCR Purification** | Beckman Coulter Life Sciences, Catalog #A63881 |
| **Illumina DNA Prep library kit** | Illumina, Catalog #20018705 (96 samples), or #20018704 (24 samples) |
| **\*IDT® for Illumina® DNA UD Index kits (plate)** | Illumina, Catalog #20027213 (Index set A), #20027214 (Index set B), #20027215 (Index set C), and #20027216 (Index set D).|
| **200 Proof Ethanol** | Decon Labs, Inc. Catalog #2716 |
| **AlumaSeal II aluminum seals** | Excel Scientific, Inc. Catalog #AF100 |
| **Clear, 8-strip PCR tubes domed caps** | LabSource, Catalog #T54-203-CS/10PK MFG# - 321-10-062 |
| **[Optional] 96‐well storage plates, round well, 300 uL ("MIDI" plate)** | Fisher Scientific, Catalog #AB‐0859 |
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

1. **PET-PCR Sample Quality Check**
Real-time PCR hands-on time 30 min / 96 samples; Cycle time 1.2 hours
Reagents: Primers, 2X ABI TaqMan buffer, DNase PCR free water


2. **PCR reaction to generate amplicons**
PCR hands-on time 30 min / 96 samples; Cycle time 2.5 hours
Reagents: 10uM Primers, HF Phusion Taq, 5X GC Buffer, 10mM dNTPs, DNase PCR free water


3. **Analysis of PCR amplicons**
PCR amplicon electrophoresis hands-on time 10 min / 8 samples; Gel running time 30 min
Reagents: Agarose, DNA loading dye, 1kb DNA ladder, 1X TBE Buffer
_If <20 samples, run all samples on the gel; If >20 samples, pick 20 samples with varying CT values and run on the gel_

4. **PCR amplicons clean up**
Hands on time 30 min / 96 samples; Total time 40+ min / 96 samples
Reagents: AMPure XP beads, fresh 70% EtOH, Nuclease-free water

5. **Tagment Genomic DNA and Tagmentation Clean-Up**
Hands on time 30 min / 96 samples; Total time 17 min / 8 samples
Reagents: BLT, TB1, TSB, TWB
_[optional]_ To assess tagmentation, run 1 uL sample on Agilent Bioanalyzer 2X and/or TapeStation 2X using High Sensitivity DNA chip

6. **Amplification of Tagmented DNA (Index PCR)** <br>
Hands on time 35 min / 96 samples; Cycle time 38 min / 96 samples
Reagents: EPM, Nuclease-free water, Index 1 and 2 primers

7. **Library PCR Clean-up**
Hands on time 30 min / 96 samples; Total time 40+ min / 96 samples
Reagents: SPB, RSB, Nuclease-free water, fresh 80% EtOH

8. **Library Pooling, Quantification, and Normalization**
Hands on time 30+ min / 96 samples; Total time 40+ min / 96 samples
Reagents: Sample Buffer, D5000 Ladder, ScreenTape; Qubit dsDNA HS Buffer and Reagent, Standard #1 and #2

9. **Library Denaturing and MiSeq Sample Loading**
Hands on time 30 min / pooled samples; Total time 30 min / pooled samples
Reagents: Resuspension Buffer, HT1, 0.2N NaOH, PhiX Control Kit v3, 200mM Tris-HCl pH7.0

10. **Analysis of NGS data**
Hands on time 5 min / 96 samples; Total time 15-25 min / 96 samples
Method: MaRS analysis pipeline

11. **Standardized SNPs reports generated**



  * * *


<a id="pet_pcr"></a>

## Sample QC ##

This step uses a real time PCR assay, [PET-PCR](https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0179178),  to assess the quality and quantity of starting DNA material. The readout includes an estimation of _all_ DNA in the sample, host and parasite.

- **NOTE:** This assay is used to identify which samples should be included and/or excluded for downstream procedures. Any sample with a CT value greater than 34, should be excluded. We recommend grouping all samples with a CT > 34 and performing separate amplicon PCRs and electrophoresis for these samples. If the electrophoresis yields positive results (e.g. positive and clear bands on gel) you may procced with downstream procedures. Samples with a CT > 34 have lower  parasite DNA and subsequent gene amplification success rates are variable.

**Consumables**

**Table 4. PET-PCR Consumables**

| Item | Quantity | Storage |
| --- | --- | --- |
| **Primers – FAM labled genus primers and HEX labeled _falciparum_ primers (see below)** | 0.25-0.5 uL per sample | 2° to 8°C |
| **TaqMan 2X Environmental buffer** | 10 uL per sample | 2° to 8°C |
| **Nuclease-free water** | 6.25 uL per sample | Room temperature |
| **Strip tubes 8X** | Up to 8 samples per strip | Room temperature |
| **Strip Optical caps 8X** | Up to 8 samples per strip | Room temperature |

**Preparation**

- All stock primers should be prepared at a 10uM concentration.<br>
- DNA samples should be stored at 4C until testing or -20C for long term storage.<br>
- Store all primer stocks at -20C for up to 1 year.<br>
- Unopened tubes of ABI TaqMan Environmental Buffer should be stored at -20C for a maximum of six months. Once thawed, store at 4C for up to six months. The reagent must be used within the expiration date provided by the manufacturer.<br>
- All samples should be tested in duplicates or triplicates in some special cases (e.g. very low density situations).<br>
- Typically, the genus _P. falciparum_ multiplex assay should be run first on all samples. All genus-positive sample are subsequently tested for _P. ovale, P. malaria_ and _P. vivax_, in order to determine the species.

**Procedure**

**Initial Set up**

- The PET-PCR reaction mix is prepared by mixing the TaqMan environmental buffer, primers, and water as shown below.
- Determine the number of reactions you need to run by multiplying the total number of samples you have to test (including your positive and negative controls) by two because every sample will be tested in duplicates.
  - For example, if you are testing 10 samples, you will multiple this by two to give you 20. Add two extra reactions to account for loss of solution during pipetting. This gives you a total number of 22 reactions. Multiply this number with the volumes below for each component to give you the total master-mix volume required for your experiment.
- In a 1.5mL tube, prepare your master-mix by multiplying the volumes shown below with the total number of reactions you need to run (e.g. 22 as described above).


**Primers and PCR Conditions**

The table below shows the Genus and _P.falciparum_ primers and PCR conditions for a multiplex reaction:  

**Table 5. Multiplexing Genus and _P. falciparum_ species specific primers**

| Master mix              | Reaction volume | x samples + 1 | Final [conc] |
|----                     |----             |----           |----          |
|  water                  | 6.25 uL         |               |              |
|  2X ABI TaqMan buffer   | 10.00 uL        |               |    1x        |
|  Genus F primer         | 0.50 uL         |               |    0.250 uM  |
|  FAM-genus R primer     | 0.50 uL         |               |    0.250 uM  |
|  P.f. F primer          | 0.50 uL         |               |    0.250 uM  |
|  HEX-P.f. R primer      | 0.25 uL         |               |    0.125 uM  |
|  **TOTAL**              | **18.0 uL**     |               |              |
|  **Add last**           |                 |               |              |
|  DNA                    | 2.0 uL          |               |              |

>Primers:

>**Genus 18sFor (5' to 3'):** 5'-GGC CTA ACA TGG CTA TGA CG-3'

>**Genus 18sRev (5' to 3'):** 5'-agg cgc ata gcg cct ggC TGC CTT CCT TAG ATG TGG TAG CT-3' (**FAM-labeled: based on the 18s rRNA gene)**

>_**P. falciparum**_ **For (5' to 3'):** 5'-ACC CCT CGC CTG GTG TTT TT-3'

>_**P. falciparum**_ **Rev (5' to 3'):** 5'-agg cgg ata ccg cct ggT CGG GCC CCA AAA ATA GGA A-3' **(HEX-labeled: based on the r364 target)**

**Thermocyclying conditions**:

For 1 cycle:

| Temperature    | Time (min)    |
|----            |----           |
| 95&deg;C       | 15:00         |      

Then, for 45 cycles:

| Temperature    | Time (min)    |
|----            |----           |
| **95&deg;C**   | **0:20**      |        
| **63&deg;C**   | **0:40**      |        
| **72&deg;C**   | **0:30**      |       

Lastly, hold:

| Temperature    | Time (min)    |
|----            |----           |
| 4&deg;C        | ∞             |

**Adding the DNA Samples**

1. Mix the prepared master-mix well by vortexing briefly.

2. Centrifuge the tubes for 5 seconds to remove any solution trapped in the cap.

3. Arrange the optically clear PCR tubes on a PCR-tube rack following the PCR sample sheet. Add 18 uL of the PET-PCR master mix prepared above to each PCR well. Loosely put on the lids of the wells filled with master mix solution.

4. Return all reagents to the freezer and refrigerator before proceeding to the next step.

5. Take the assembled plate containing the tubes with PCR master mix solution to the PCR template area.

6. Add 2 uL of the unknown DNA samples to the wells with the master-mix according to the sample sheet. Cap the well tightly after adding the sample. The total volume of PCR reaction is 20.0 uL after addition of the template.

7. Add positive control DNA to each positive control well with master-mix. Cap the wells after each positive control is added.

8. Add 2.0 uL of DNase-free water to the wells designated as the no-template control (NTC) and close that well tightly.

9. Make sure each sample has been added to the correct well and that all wells are tightly capped.

10. Briefly centrifuge your strip tubes to remove any solution trapped on the walls of the wells.

11. Make sure there are no bubbles in the well.

- **NOTE:** _The recommended minimum amount of template DNA is 2.0 uL. This can be adjusted appropriately depending on the sample parasitemia._

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
Lucchi, N.W., et al., _Molecular diagnosis of malaria by photo-induced electron transfer fluorogenic primers: PET-PCR._ PLoS One, 2013. 8 (2): p. e56677.


  * * *


## Gene Enrichment & QC ##

#### [I: Gene PCR Amplification](#gene_amp) ####
#### [II: QC by Electrophoresis](#electrophoresis) ####

### Recommended Positive Controls and Expected SNPs ###

We highly recommend using at least **three positive controls** with known SNP profiles. These will be first analyzed to confirm the known SNPs and ensure the rest of the sequencing run was successful.

We routinely use the following controls:

| Control strain | _crt_         | _mdr1_        | _dhfr_        | _dhps_    | _k13_       |
| ---            | ---           | ---           | ---           | ---       | ---         |
| 7G8            | **S**VMN**T** | NED**FCDY**   | C**I**C**N**I | S**G**KAA | _wild type_ |
| DD2            | CV**IET**     | **Y/F**EDFCDY | C**IRN**I     | S**G**KAA | _wild type_ |
| HB3           
