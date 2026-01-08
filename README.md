## Table of contents
1. [About this Repository](#About-this-Repository)
2. [Overview](#Overview)
3. [Step_1–2_LDSC_estimation](#Step_1_2_LDSC_estimation)
4. [Step_3_PHBC_modelling](#Step_3_PHBC_modelling)
5. [PRS_CS](#PRS_CS)
4. [UKB_PheWAS](#UKB_PheWAS)
5. [All_of_Us_PheWAS](#All_of_Us_PheWAS)
6. [Nigthshift_stratified](#Nigthshift_stratified)

## About this Repository
This repository accompanies the manuscript __"Insights into Pathophysiological Pathways in ME/CFS Through Genetic Correlation and Mendelian Randomization"__
<br/><br/> <br/>
<p align="center">
<img src="/img/Fig_1_plusPHBC_final.png" alt="Figure 1" width="600"/>
<br/><br/>



> **_Abstract:_**  Myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS) and post-acute infection syndromes (PAIS) are multisystem disorders in which immune, vascular, neuroinflammatory, and metabolic disturbances are frequently described, yet their causal roles remain unclear. Leveraging genome-wide summary statistics from the DecodeME study 15579 cases), we used genetic correlation, shared pleiotropic heritability, and Mendelian randomization to test whether these proposed biological domains have a genetic basis in ME/CFS. Across 22 auxiliary traits, representing five broad mechanistic domains, we observed significant positive genetic correlations for traits of cellular energetics, neurovascular regulation, and epithelial barrier–microbiome function. Migraine and irritable bowel syndrome, proxies for neurovascular dysregulation and barrier–microbiome disturbance, showed the strongest genetic overlap with ME/CFS, each contributing substantially to shared pleiotropy. In contrast, platelet/immunothrombotic and inflammatory traits showed smaller but still measurable genetic correlations and pleiotropic contributions. Energetics-related traits, including type 2 diabetes and lactate, displayed consistent genetic correlation but comparably low shared pleiotropy, suggesting that metabolic inflexibility may act through broader physiological networks. Mendelian randomization analyses identified three biomarkers with evidence for a causal contribution to ME/CFS: higher mitochondrial DNA copy number was protective, whereas increased glycoprotein acetyls (systemic inflammation) and mean platelet volume (platelet activation) increased disease risk. Together, these results show that ME/CFS susceptibility is shaped by interacting pathways involving barrier–microbiome dysfunction, neurovascular instability, chronic inflammation, platelet activation, and impaired cellular energetics. This genetically anchored framework highlights candidate mechanisms for biomarker and therapeutic development in ME/CFS and PAIS.
<p>
<br/>


This repository contains a reproducible, HPC-friendly implementation of the   **Pleiotropic Heritability Block Correction (PHBC)** framework and was used to generate **Figure 1B**  

__The implementation is based on the methodology and core code described in:__
>* [__pleioh2g preprint:__](https://doi.org/10.1101/2025.06.10.25329261) https://doi.org/10.1101/2025.06.10.25329261  
>* [__Original code:__](https://github.com/cran/pleioh2g) https://github.com/cran/pleioh2g  


## Overview
PHBC extends LD Score Regression by partitioning SNP-based heritability into trait-specific and shared pleiotropic components across multiple traits.This repository provides a **modular two-step pipeline** that applies PHBC at scale, enabling leave-one-out and multi-trait analyses in distributed computing environments.  

Because the original *pleioh2g* R package was difficult to install and the wrapper functions did not run reliably in our environment, we implemented a lightweight pipeline that directly uses the underlying functions while adding:
>* Robust job submission on HPC/SLURM systems  
>* Explicit control over LDSC and jackknife steps  
>* Reproducible software [__deployment via Docker__](https://hub.docker.com/r/mwielsch/pleioh2g)    



## Step_1_2_LDSC_estimation 

[__Script__](code/PHBC_step1_2_ldsc.R)  

>* **Purpose:**  
>>* Computes SNP-based heritability and pairwise genetic correlations between a target trait (e.g. ME/CFS) and a panel of auxiliary traits using LDSC.

>* **Inputs:**
>>* Pre-munged GWAS summary statistics (`*.sumstats.gz`)
>>* LD score reference files
>>* Sample and population prevalence estimates (for liability-scale correction)

>* **Outputs:**
>>* LDSC heritability and genetic correlation estimates
>>* Jackknife-based standard errors
>>* Intermediate `.rds` files used by PHBC
>>* Summary tables used for downstream analyses and figures

This step produces the genetic correlation estimates shown in **Figure 1A** of the manuscript.


## Step_3_PHBC_modelling 

[__Script__](code/PHBC_step3_pleioh2g_biofunctions.R)  

>* **Purpose:**  
>>* Applies the PHBC model to decompose trait heritability into:
>>* Trait-specific components
>>* Shared pleiotropic components across traits

>* **Analyses performed:** 
>>* Single-trait PHBC models (pairwise with the target trait)
>>* Leave-one-out models excluding each auxiliary trait in turn
>>* A full multi-trait PHBC model including all traits

>* **Outputs:**
>>* Estimates of shared pleiotropic heritability
>>* Leave-one-out differences used to infer domain-specific pleiotropic contributions

These results are summarized in **Figure 1B** of the manuscript.





