#!/bin/bash
#SBATCH -J ME_CFS_step1
#SBATCH -p zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --tasks-per-node=5
#SBATCH --mem=95G


module load singularity

singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
  Rscript PHBC_step1_2.R ME-CFS

singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
 Rscript PHBC_step1_2.R ME-CFS-female

singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
 Rscript PHBC_step1_2.R ME-CFS-male

singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
 Rscript PHBC_step1_2.R ME-CFS-infect

singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
 Rscript PHBC_step1_2.R ME-CFS-nonInfect


singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
 Rscript PHBC_step1_2.R LC

#### ----------------------------------------------------
### submit second part of pipeline
###   ---------------------------------------------

#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R ME-CFS

#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R ME-CFS-female

#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R ME-CFS-male

#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R ME-CFS-infect


#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R ME-CFS-nonInfect

#singularity exec -B /gpfs:/gpfs pleioh2g_v1.0.sif \
#  Rscript PHBC_step3_pleioh2g_biofunctions.R LC



