#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tibble)
library(fs)
library(pleioh2g)

## ---- read target from command line (default: ME-CFS) ----
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "ME-CFS"
message("Running LDSC for target: ", target)

## ---- paths & traits ----
base_dir   <- "/gpfs/data/fs72678/mwielsch3/ME_CFS/MR_energy/PHBC"
sumdir     <- file.path(base_dir, "input")
result_dir <- file.path(base_dir, "result")
tmp_dir    <- file.path(base_dir, "tmp")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir,    showWarnings = FALSE, recursive = TRUE)

ld_path  <- "/gpfs/data/fs72678/mwielsch3/ME_CFS/LDscoreREG/eur_w_ld_chr"
wld_path <- "/gpfs/data/fs72678/mwielsch3/ME_CFS/LDscoreREG/eur_w_ld_chr"
hmp3_local <- "/gpfs/data/fs72678/mwielsch3/ME_CFS/LDscoreREG/w_hm3.snplist"
hmp3 <- if (file_exists(hmp3_local)) hmp3_local else
  fs::path(fs::path_package("extdata/w_hm3.snplist", package = "pleioh2g"))

## ---- lead trait & auxiliary panel ----
aux_traits <- c(
  "WMH-volume",
  "Alzheimers",
  "Mito-DNA-CN",
  "Varicose-veins",
  "rbc-mono",
  "rbc-meanPlatelet-volume",
  "VWF",
  "MS",
  "Crohns",
  "nmr-Glycoprotein-Acetyls",
  "Venous-thromboembolism",
  "Ulcerative-colitis",
  "T2D",
  "nmr-Lactate",
  "SLE",
  "ortho-hypotens",
  "Cortisol-level-1",
  "Insomina",
  "Migrane",
  "IBS"
)

phenotype <- c(target, aux_traits)
G <- which(phenotype == target)

## used for filenames; replace non-alnum with "_"
target_tag <- gsub("[^A-Za-z0-9]+", "_", target)

## ---- prevalence lookup ----
## sample_prev = case / (case + control) in the GWAS sample
## population_prev = population prevalence on liability scale
## (fill in / adjust values as needed)

sample_prev_lookup <- c(
  "ME-CFS"              = 0.05655,
  "ME-CFS-female"       = 0.05655,  
  "ME-CFS-male"         = 0.05655,  
  "ME-CFS-infect"       = 0.05655,  
  "ME-CFS-nonInfect"    = 0.05655,  
  "LC"                  = 0.00588, 
  "Alzheimers"          = 0.08019,
  "Varicose-veins"      = 0.08223,
  "MS"                  = 0.40957,
  "Crohns"              = 0.06513,
  "Ulcerative-colitis"  = 0.28195,
  "T2D"                 = 0.08253,
  "SLE"                 = 0.02759,
  "ortho-hypotens"      = 0.00346,
  "Insomina"            = 0.26576,
  "Migrane"             = 0.10,     # approx
  "IBS"                 = 0.10974
)

population_prev_lookup <- c(
  "ME-CFS"              = 0.005,
  "ME-CFS-female"       = 0.005,    
  "ME-CFS-male"         = 0.005,
  "ME-CFS-infect"       = 0.005,    
  "ME-CFS-nonInfect"    = 0.005,
  "LC"                  = 0.027, 
  "Alzheimers"          = 0.10,
  "Varicose-veins"      = 0.157,
  "MS"                  = 0.00375,
  "Crohns"              = 0.003,
  "Ulcerative-colitis"  = 0.0023,
  "T2D"                 = 0.10,
  "SLE"                 = 0.003,
  "ortho-hypotens"      = 0.22,
  "Insomina"            = 0.124,
  "Migrane"             = 0.12,
  "IBS"                 = 0.092
)

## helper: get prevalences with warning if missing
get_prev_vector <- function(phenos, lookup, what = "sample") {
  v <- unname(lookup[phenos])
  missing <- is.na(v) & !(phenos %in% names(lookup))
  if (any(missing)) {
    warning(
      "Missing ", what, " prevalence for: ",
      paste(phenos[missing], collapse = ", "),
      ". Using NA for these."
    )
  }
  v
}

sample_prev     <- get_prev_vector(phenotype, sample_prev_lookup, what = "sample")
population_prev <- get_prev_vector(phenotype, population_prev_lookup, what = "population")

## ---- read sumstats ----
read_sumstats_local <- function(path, name) {
  message("Reading ", name, " from ", path)
  df <- read_tsv(path, col_types = cols())
  df <- dplyr::filter(df, !is.na(Z))
  as_tibble(df)
}

file_for <- function(tr) fs::path(sumdir, sprintf("munge_1_%s.sumstats.gz", tr))
files <- setNames(vapply(phenotype, file_for, character(1)), phenotype)

if (!all(file_exists(files))) {
  stop("Missing sumstats files for: ",
       paste(names(files)[!file_exists(files)], collapse = ", "))
}

munged_sumstats <- lapply(phenotype, function(tr) {
  read_sumstats_local(files[[tr]], tr)
})
names(munged_sumstats) <- phenotype

## ============================================================
## Step 1: LDSC (no jackknife)
## ============================================================
message("Step 1: Cal_rg_h2g_alltraits for target ", target, "...")

ldsc_res <- tryCatch(
  Cal_rg_h2g_alltraits(
    phenotype       = phenotype,
    munged_sumstats = munged_sumstats,
    ld_path         = ld_path,
    wld_path        = wld_path,
    sample_prev     = sample_prev,
    population_prev = population_prev
  ),
  error = function(e) {
    message("Cal_rg_h2g_alltraits failed: ", conditionMessage(e))
    stop(e)
  }
)

ldsc_main_path <- file.path(tmp_dir, paste0("ldsc_", target_tag, "_main.rds"))
saveRDS(ldsc_res, ldsc_main_path)

## quick summary vs target
rg_vec  <- as.numeric(ldsc_res$rg[target, ])
rgz_vec <- as.numeric(ldsc_res$rgz[target, ])
rg_se   <- ifelse(rgz_vec != 0, rg_vec / rgz_vec, NA_real_)
rg_p    <- 2 * pnorm(-abs(rgz_vec))

h2_obs    <- as.numeric(ldsc_res$h2[1, ])
h2_obs_se <- ifelse(ldsc_res$h2Z[1, ] != 0, h2_obs / ldsc_res$h2Z[1, ], NA_real_)
h2_lia    <- if (!is.null(ldsc_res$liah2)) as.numeric(ldsc_res$liah2[1, ]) else rep(NA_real_, length(phenotype))

ldsc_tbl <- tibble(
  p1          = phenotype,
  target_trait= target,
  rg          = rg_vec,
  rg_se       = rg_se,
  rg_z        = rgz_vec,
  rg_p        = rg_p,
  h2_obs      = h2_obs,
  h2_obs_se   = h2_obs_se,
  h2_lia      = h2_lia,
  sample_prev = sample_prev,
  pop_prev    = population_prev
) %>% filter(p1 != target)

ldsc_csv_path <- file.path(result_dir, paste0("ldsc_", target_tag, "_summary.csv"))
write_csv(ldsc_tbl, ldsc_csv_path)

message("Step 1 done. Saved:")
message("  ", ldsc_main_path)
message("  ", ldsc_csv_path)

## ============================================================
## Step 2: LDSC jackknife
## ============================================================
message("Step 2: Cal_rg_h2g_jk_alltraits (jackknife) for target ", target, "...")

n_block <- 200  ## adjust if needed

ldsc_jk <- tryCatch(
  Cal_rg_h2g_jk_alltraits(
    n_block         = n_block,
    hmp3            = hmp3,
    phenotype       = phenotype,
    munged_sumstats = munged_sumstats,
    ld_path         = ld_path,
    wld_path        = wld_path,
    sample_prev     = sample_prev,
    population_prev = population_prev
  ),
  error = function(e) {
    message("Cal_rg_h2g_jk_alltraits failed: ", conditionMessage(e))
    stop(e)
  }
)

ldsc_jk_path <- file.path(tmp_dir, paste0("ldsc_", target_tag, "_jk.rds"))
saveRDS(ldsc_jk, ldsc_jk_path)

message("Step 2 done. Saved:")
message("  ", ldsc_jk_path)
















