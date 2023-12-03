#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --job-name=rnb
#SBATCH --mail-type=ALL
#SBATCH --output ./log.txt
#SBATCH --mail-user=ccthomas10@gmail.com

# Load Packages and Dependencies
ml palma/2020b
ml GCC/10.2.0
ml OpenMPI/4.0.5
ml DecompPipeline/1.0.3-R-4.0.3

# Call Decomp Pipeline
Rscript run_rnb.R
