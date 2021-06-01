#!/bin/sh

#SBATCH -J SLiM_Simulation
#SBATCH -t 128:00:00
#SBATCH -p apophis
#SBATCH -o sim.out
#SBATCH -e sim.err
#SBATCH -n 1

for i in {1..100}; do mkdir rep_$i; cp one_point_five_subs_0.2.tree rep_$i; cp run_simulation.sh rep_$i; cd rep_$i; sbatch run_simulation.sh; cd ..;  done 
