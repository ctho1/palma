#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --partition normal
#SBATCH --time=8:00:00
#SBATCH --mem=20G
#SBATCH --job-name=downloader
#SBATCH --mail-type=ALL
#SBATCH --output log.txt
#SBATCH --mail-user=christian.thomas@ukmuenster.de

cat download-links.txt | parallel --gnu "wget --content-disposition {}"