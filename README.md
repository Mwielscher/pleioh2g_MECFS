## Table of contents
1. [About this Repository](#About-this-Repository)
2. [Overview](#Overview)
3. [MultiVar_GWAS](#MultiVar_GWAS)
4. [Mendelian_Randomisation](#Mendelian_Randomisation)
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
<p>
Because the original *pleioh2g* R package was difficult to install and the wrapper functions did not run reliably in our environment, we implemented a lightweight pipeline that directly uses the underlying functions while adding:
>* Robust job submission on HPC/SLURM systems  
>* Explicit control over LDSC and jackknife steps  
>* Reproducible software [__deployment via Docker__](https://hub.docker.com/r/mwielsch/pleioh2g)    



## MultiVar_GWAS  
We incorporated 21 sets of summary statistics, totaling 3.5 million SNPs, into our structural multivariable regression model. All models were based on samples of European ancestry or trans-ethnic meta-analyses, with European linkage disequilibrium maps as references.
To mitigate type 1 error inflation, we used the built-in genomic control function of the Genomic SEM package. Additionally, we employed the "smooth_check" option to exclude SNPs requiring excessive smoothing, enhancing the reliability of the GWAS results. The "fix_measurement" option was also applied to ensure stable factor loadings, preventing implausible estimates and Heywood cases.
__GCTA:__ We used the UK Biobank genotyping array data as our reference sample (see UKB genetic data QC section). A genome-wide significance threshold of 5x10⁻⁸ was applied. We analyzed a 10 MB window around each SNP, assuming that SNPs outside this window are in complete linkage equilibrium. Additionally, we applied a collinearity threshold of 0.9 to exclude highly correlated SNPs.


## Mendelian_Randomisation
Scripts include: 
__Proxy search__ - we used proxies in high linkage disequilibrium (LD) with the original instruments, applying an LD threshold of 0.8 in the European sample of the 1000 Genomes Project Phase 3. Proxy SNPs were retrieved using the Ensembl REST API. We then selected the SNP with the lowest P-value from our summary statistics.
__Pleiotropy evaluation;__ we also calculate heterogeneity indicators specific to each outcome, including Cochran's Q test results, which assess variability in causal estimates and potential horizontal pleiotropy. We also evaluated the Egger intercept to flag directional pleiotropy in the instruments. Scripts also return Steiger P values for directionality  
__Sensitivity analysis:__ radial IVW method, designed to handle outliers within the instrument set by re-weighting the contribution of each variant based on its residual, making it more robust to heterogeneity among genetic instruments. Additionally, we included MR-RAPS to handle cases where some genetic instruments are weak, mitigating weak instrument bias.
Finally, we reported results from the MR-PRESSO method, including the MR-PRESSO Global P-value, which tests for horizontal pleiotropy by simulating a null distribution of residual sum of squares and comparing it to the observed values. This identifies whether the observed residual sum of squares falls within the expected variance or suggests horizontal pleiotropy. We also provided the MR-PRESSO Distortion P-value, which compares the causal estimate before and after outlier correction, indicating the influence of outliers on the causal estimate.


## PRS_CS
We calculated scores using the PRS-CS approach, which employs a Bayesian framework that incorporates linkage disequilibrium (LD) information from large reference panels. This method adaptively shrinks SNP effect sizes, allowing for the flexible handling of both large and small genetic effects.  
__MTAG__ We ran MTAG using its default settings, thus obtained the same number of summary statistics as the input traits, with each set containing updated meta-analyzed effect estimates based on the genetic correlations among the included traits. We then focused on the MTAG results where effect estimates were driven by the sleep factor, rather than the psychiatric or behavioral traits. These adjusted factors were subsequently used in PheWAS analysis and polygenic score (PGS) analysis, conducted separately in night shift and day shift workers.  

## UKB_PheWAS
The following data fields were used: 20002 (Self-reported non-cancer diagnosis), 20008 (Year of self-reported diagnosis), 20004 (Self-reported operation code), 40006 (Cancer diagnosis from cancer registry), and 40005 (Date of cancer diagnosis). Additionally, we retrieved record-level data from 42040 (GP clinical event records), 42039 (GP prescription records), 41259 (Hospital Episode Statistics [HES] inpatient main dataset), 41234 (HES inpatient diagnosis), 41149 (HES inpatient operations), 40023 (Cause of death), and 40023 (Date of death). Using DeepPheWAS, we mapped these data to phecodes, created quantitative phenotypes, and generated combined phenotypes. This resulted in a dataset containing 1,234 variables, each with at least 20 entries in the analysed UKB subpopulation 
