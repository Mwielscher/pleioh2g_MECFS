#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tibble)
library(fs)
library(pleioh2g)

## ============================================================
## 1. Target from command line
## ============================================================
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "ME-CFS"
message("Step 3 (per trait): running pleioh2g for target: ", target)

target_tag <- gsub("[^A-Za-z0-9]+", "_", target)

## ============================================================
## 2. Paths
## ============================================================
base_dir   <- "/gpfs/data/fs72678/mwielsch3/ME_CFS/MR_energy/PHBC"
result_dir <- file.path(base_dir, "result")
tmp_dir    <- file.path(base_dir, "tmp")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## 3. Map traits -> category  (for later plotting)
## ============================================================

trait_category <- c(
  # Energy / Cellular energetics
  "T2D"                   = "Energy_Cellular_energetics",
  "nmr-Total-Triglycerides" = "Energy_Cellular_energetics",
  "nmr-Lactate"           = "Energy_Cellular_energetics",
  "Mito-DNA-CN"           = "Energy_Cellular_energetics",

  # Immunothrombosis / Platelets
  "rbc-meanPlatelet-volume" = "Immunothrombosis_Platelets",
  "VWF"                      = "Immunothrombosis_Platelets",
  "Venous-thromboembolism"   = "Immunothrombosis_Platelets",

  # Perfusion / Vascular
  "Varicose-veins"        = "Perfusion_Vascular",
  "ortho-hypotens"        = "Perfusion_Vascular",
  "Migrane"               = "Perfusion_Vascular",
  "WMH-volume"            = "Perfusion_Vascular",

  # Barrier / Microbiome–epithelium
  "IBS"                   = "Barrier_Microbiome_epithelium",
  "Crohns"                = "Barrier_Microbiome_epithelium",
  "Ulcerative-colitis"    = "Barrier_Microbiome_epithelium",

  # Neuroinflammation & Adaptive skew
  "Alzheimers"            = "Neuroinflammation_Adaptive_skew",
  "MS"                    = "Neuroinflammation_Adaptive_skew",
  "Cortisol-level-1"      = "Neuroinflammation_Adaptive_skew",
  "Insomina"              = "Neuroinflammation_Adaptive_skew",
  "SLE"                   = "Neuroinflammation_Adaptive_skew",
  "RheumatoidArthritis"   = "Neuroinflammation_Adaptive_skew",

  # Inflammatory
  "nmr-Glycoprotein-Acetyls" = "Inflammatory",
  "rbc-mono"                 = "Inflammatory"
)

all_mapped_traits <- names(trait_category)

## ============================================================
## 4. Read LDSC outputs (from step1_2)
## ============================================================
ldsc_res_path <- file.path(tmp_dir, paste0("ldsc_", target_tag, "_main.rds"))
ldsc_jk_path  <- file.path(tmp_dir, paste0("ldsc_", target_tag, "_jk.rds"))

if (!file_exists(ldsc_res_path) || !file_exists(ldsc_jk_path)) {
  stop(
    "Missing LDSC RDS files for target ", target,
    ". Did PHBC_step1_2_ldsc.R run successfully for this target?"
  )
}

ldsc_res <- readRDS(ldsc_res_path)
ldsc_jk  <- readRDS(ldsc_jk_path)

ldsc_pheno <- colnames(ldsc_res$rg)

## ============================================================
## 5. Define available traits
## ============================================================
missing_global <- setdiff(all_mapped_traits, ldsc_pheno)

if (length(missing_global) > 0) {
  warning(
    "The following mapped traits are not in LDSC results and will be ignored: ",
    paste(missing_global, collapse = ", ")
  )
}

available_traits <- intersect(all_mapped_traits, ldsc_pheno)

message("Available auxiliary traits in LDSC for target ", target, ":")
message("  ", paste(available_traits, collapse = ", "))

## ============================================================
## 6. Helper function: run pleioh2g for a given aux set
## ============================================================
run_pleio_for_set <- function(aux_traits,
                              unit_name,    # here: trait name
                              mode = c("only", "leave_one_out"),
                              target,
                              target_tag,
                              ldsc_res,
                              ldsc_jk,
                              result_dir,
                              ldsc_pheno,
                              trait_category) {

  mode <- match.arg(mode)
  trait_tag <- gsub("[^A-Za-z0-9]+", "_", unit_name)
  mode_tag  <- mode

  message("-------------------------------------------------")
  message("Trait: ", unit_name, " | mode: ", mode)

  ## filter aux_traits to those present in LDSC
  existing_aux <- intersect(aux_traits, ldsc_pheno)
  removed_aux  <- setdiff(aux_traits, existing_aux)

  if (length(removed_aux) > 0) {
    message("  [INFO] Removing missing traits: ", paste(removed_aux, collapse = ", "))
  }

  if (length(existing_aux) == 0L) {
    warning("  No usable aux traits left for ", unit_name, " (mode=", mode, "). Skipping.")
    return(NULL)
  }

  phenotype <- c(target, existing_aux)
  G <- which(phenotype == target)
  if (length(G) != 1) stop("Target phenotype missing after filtering.")

  ## h2 and rg extraction
  if (!is.null(ldsc_res$liah2)) {
    h2_vector <- as.numeric(ldsc_res$liah2[1, phenotype])
    h2_matrix <- ldsc_jk$liah2array[, phenotype]
  } else {
    h2_vector <- as.numeric(ldsc_res$h2[1, phenotype])
    h2_matrix <- ldsc_jk$h2array[, phenotype]
  }

  Results_full_rg       <- ldsc_res$rg[phenotype, phenotype]
  Results_full_rg_array <- ldsc_jk$rgarray[phenotype, phenotype, ]

  ## ---------- pleioh2g pre-correction ----------
  message("  pleiotropyh2_nocor_computing_single ...")
  pleio_precorr <- tryCatch(
    pleiotropyh2_nocor_computing_single(
      G                     = G,
      phenotype             = phenotype,
      h2_vector             = h2_vector,
      h2_vector_mat         = h2_matrix,
      Results_full_rg       = Results_full_rg,
      Results_full_rg_array = Results_full_rg_array
    ),
    error = function(e) {
      message("  [WARN] precorr failed: ", conditionMessage(e))
      NULL
    }
  )

  ## ---------- pleioh2g bias-corrected ----------
  message("  pleiotropyh2_cor_computing_single_prune ...")
  sample_rep <- 500

  pleio_corr <- tryCatch(
    pleiotropyh2_cor_computing_single_prune(
      G                     = G,
      phenotype             = phenotype,
      h2_vector             = h2_vector,
      h2_vector_mat         = h2_matrix,
      Results_full_rg       = Results_full_rg,
      Results_full_rg_array = Results_full_rg_array,
      sample_rep            = sample_rep
    ),
    error = function(e) {
      message("  [WARN] corrected run failed: ", conditionMessage(e))
      NULL
    }
  )

  ## ---------- save raw objects ----------
  prefix <- file.path(
    result_dir,
    paste0("pleioh2g_", target_tag, "_", trait_tag, "_", mode_tag)
  )

  saveRDS(pleio_precorr, paste0(prefix, "_precorr.rds"))
  saveRDS(pleio_corr,    paste0(prefix, "_corr.rds"))

  ## ---------- summary table ----------
  extract <- function(obj, field) {
    if (is.null(obj)) return(NA_real_)
    v <- obj[[field]]
    if (is.null(v)) return(NA_real_)
    as.numeric(v)
  }

  category <- trait_category[unit_name]
  if (is.na(category)) category <- NA_character_

  summary_tbl <- tibble(
    target                        = target,
    aux_trait                     = unit_name,
    category                      = category,
    mode                          = mode_tag,
    auxiliaries_input             = paste(aux_traits, collapse = ";"),
    auxiliaries_used              = paste(existing_aux, collapse = ";"),
    h2_target                     = extract(pleio_precorr, "target_disease_h2_est"),
    h2_target_se                  = extract(pleio_precorr, "target_disease_h2_se"),
    h2_pleiotropic_precorr        = extract(pleio_precorr, "h2pleio_uncorr"),
    h2_pleiotropic_precorr_se     = extract(pleio_precorr, "h2pleio_uncorr_se"),
    pct_h2_pleiotropic_precorr    = extract(pleio_precorr, "percentage_h2pleio_uncorr"),
    pct_h2_pleiotropic_precorr_se = extract(pleio_precorr, "percentage_h2pleio_uncorr_se"),
    h2_pleiotropic_corr           = extract(pleio_corr, "h2pleio_corr"),
    h2_pleiotropic_corr_se        = extract(pleio_corr, "h2pleio_corr_se"),
    pct_h2_pleiotropic_corr       = extract(pleio_corr, "percentage_h2pleio_corr"),
    pct_h2_pleiotropic_corr_se    = extract(pleio_corr, "percentage_h2pleio_corr_se")
  )

  summary_path <- paste0(prefix, "_summary.csv")
  write_csv(summary_tbl, summary_path)

  message("  Saved: ", summary_path)

  invisible(summary_tbl)
}

## ============================================================
## 7. Loop over *traits*: ONLY + LEAVE-ONE-OUT
## ============================================================
all_summaries <- list()

for (tr in available_traits) {

  ## ===== ONLY this trait as auxiliary =====
  s_only <- run_pleio_for_set(
    aux_traits     = tr,
    unit_name      = tr,
    mode           = "only",
    target         = target,
    target_tag     = target_tag,
    ldsc_res       = ldsc_res,
    ldsc_jk        = ldsc_jk,
    result_dir     = result_dir,
    ldsc_pheno     = ldsc_pheno,
    trait_category = trait_category
  )

  ## ===== LEAVE-THIS-TRAIT-OUT: all others as auxiliary =====
  traits_leave <- setdiff(available_traits, tr)

  s_leave <- run_pleio_for_set(
    aux_traits     = traits_leave,
    unit_name      = tr,
    mode           = "leave_one_out",
    target         = target,
    target_tag     = target_tag,
    ldsc_res       = ldsc_res,
    ldsc_jk        = ldsc_jk,
    result_dir     = result_dir,
    ldsc_pheno     = ldsc_pheno,
    trait_category = trait_category
  )

  all_summaries[[paste0(tr, "_only")]]         <- s_only
  all_summaries[[paste0(tr, "_leave_one_out")]] <- s_leave
}

message("Step 3 (per trait) completed for target: ", target)









