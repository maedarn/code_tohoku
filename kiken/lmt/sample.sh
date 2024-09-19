#!/bin/sh
#------ pjsub option --------#
#SBATCH -J Zm4V100D1B1A15
#SBATCH -N 5
#SBATCH -n 512
#SBATCH -c 1
#SBATCH -p M
#------- Program execution -------#

cd /sc/home/ryunosuke.maeda/data/Zm4V100D1B1A15/
srun  ./Zm4V100D1B1A15.out > log.txt
