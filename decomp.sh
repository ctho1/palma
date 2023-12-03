#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --partition normal
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --job-name=decomp
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=ccthomas10@gmail.com

# Load Packages and Dependencies
ml palma/2020b GCC/10.2.0 OpenMPI/4.0.5 DecompPipeline/1.0.3-R-4.0.3

# Call Decomp Pipeline
Rscript decomp.R
