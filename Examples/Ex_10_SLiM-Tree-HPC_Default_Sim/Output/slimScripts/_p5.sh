#!/bin/sh

#SBATCH -J SLiM_Simulation_p5
#SBATCH -t 12:00:00
#SBATCH -p apophis
#SBATCH -o p5.out
#SBATCH -e p5.err
#SBATCH -n 1

slim /home/erin.brintnell/Summer_2021/SLiMTreeWork/slim-tree/Examples/Ex_10_SLiM-Tree-HPC_Default_Sim/Output/slimScripts/_p5.slim