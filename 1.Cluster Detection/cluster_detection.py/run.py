#!/bin/bash
#SBATCH -p short
#SBATCH --job-name="cluster_detection"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=5:00:00
#SBATCH --output=/mnt/bioinfnas2/immunocomp/clarav/analysis/hotspots/all_chr/outputs/out.%x.%j
#SBATCH --error=/mnt/bioinfnas2/immunocomp/clarav/analysis/hotspots/all_chr/errors/err.%x.%j
#SBATCH --no-requeue

echo "*"
echo "** Running cluster_detection.py for all chromosomes"
echo "*"
date
cat run_cluster_detection.sh
echo "*"
#conda activate
source /home/claravalles@vhio.org/.bashrc
source activate basics

# Command line to run the Python script for all chromosomes
python /mnt/bioinfnas2/immunocomp/clarav/scripts/Dataframes/cluster_detection.py /mnt/bioinfnas2/immunocomp/anadueso/blueprint-immuno/ANALYSIS/presentation_hotspot/ /mnt/bioinfnas2/immunocomp/clarav/analysis/hotspots/all_chr/Data/clusters_all_chr_001.tsv 001

# Compress the output file
gzip /mnt/bioinfnas2/immunocomp/clarav/analysis/hotspots/all_chr/Data/clusters_all_chr_001.tsv

echo "Finish successfully"
date
