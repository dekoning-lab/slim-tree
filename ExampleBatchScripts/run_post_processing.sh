#!/bin/sh

#SBATCH -J post_process
#SBATCH -t 128:00:00
#SBATCH -p apophis
#SBATCH -o post_process.out
#SBATCH -e post_process.err
#SBATCH -n 1

module load python/anaconda-3.6-5.1.0

for i in {1..100}; do echo $i; python3 ../../PostProcessing/process_data.py -i rep_$i/one_point_five_subs_0_nuc.fasta; cd rep_$i; 
raxmlHPC-SSE3 -s one_point_five_subs_0_nuc_consensus.fasta -f e -t ../../../PostProcessing/modelling_tree.tree -p 123456 -m GTRGAMMA -n raxml_output.tree -o p9; cd ..; done
