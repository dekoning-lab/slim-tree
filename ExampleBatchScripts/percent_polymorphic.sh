#!/bin/sh

#SBATCH -J SLiM_Simulation
#SBATCH -t 128:00:00
#SBATCH -p apophis
#SBATCH -o poly.out
#SBATCH -e poly.err
#SBATCH -n 1

module load R/3.6.2

for i in {1..100}; do cd rep_$i;
echo $i;
tail -n +2 one_point_five_subs_0_nuc_consensus_p5.tsv >  one_point_five_subs_0_nuc_consensus_p5_2.tsv;
tail -n +2 one_point_five_subs_0_nuc_consensus_p6.tsv >  one_point_five_subs_0_nuc_consensus_p6_2.tsv;
tail -n +2 one_point_five_subs_0_nuc_consensus_p7.tsv >  one_point_five_subs_0_nuc_consensus_p7_2.tsv;
tail -n +2 one_point_five_subs_0_nuc_consensus_p8.tsv >  one_point_five_subs_0_nuc_consensus_p8_2.tsv;
tail -n +2 one_point_five_subs_0_nuc_consensus_p9.tsv >  one_point_five_subs_0_nuc_consensus_p9_2.tsv;
Rscript ../../../PostProcessing/countFixations.R;
cat polymorphic_percentage.csv >> ../percent_polymorphic.csv;
cd ../;
done;
