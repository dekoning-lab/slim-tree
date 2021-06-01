#!/bin/sh

#SBATCH -J SLiM_Simulation
#SBATCH -t 128:00:00
#SBATCH -p apophis
#SBATCH -o sim.out
#SBATCH -e sim.err
#SBATCH -n 1

module load python/anaconda-3.6-5.1.0

python3 ../../../SLiMTreeHPC/SLiMTree.py -i one_point_five_subs_0.2.tree -v 5e-4 -n 100 -g 5000 -b 10 -r 2.5e-8 -k 100
