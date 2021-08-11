#!/bin/sh

#SBATCH -J SLiM_Simulation_p3
#SBATCH -t 12:00:00
#SBATCH -p apophis
#SBATCH -o p3.out
#SBATCH -e p3.err
#SBATCH -n 1

slim /home/erin.brintnell/Summer_2021/SLiMTreeWork/slim-tree/Examples/Ex_10_SLiM-Tree-HPC_Default_Sim/Output/slimScripts/_p3.slim