
#Read in the consensus data file output by the python process_data.py script
one_point_five_subs_0_nuc_consensus_p5 <- read.table("one_point_five_subs_0_nuc_consensus_p5.tsv", quote="\"", comment.char="")
one_point_five_subs_0_nuc_consensus_p6 <- read.table("one_point_five_subs_0_nuc_consensus_p6.tsv", quote="\"", comment.char="")
one_point_five_subs_0_nuc_consensus_p7 <- read.table("one_point_five_subs_0_nuc_consensus_p7.tsv", quote="\"", comment.char="")
one_point_five_subs_0_nuc_consensus_p8 <- read.table("one_point_five_subs_0_nuc_consensus_p8.tsv", quote="\"", comment.char="")
one_point_five_subs_0_nuc_consensus_p9 <- read.table("one_point_five_subs_0_nuc_consensus_p9.tsv", quote="\"", comment.char="")

#Combine into a single variable
site_files = list(one_point_five_subs_0_nuc_consensus_p5, one_point_five_subs_0_nuc_consensus_p6, 
                  one_point_five_subs_0_nuc_consensus_p7, one_point_five_subs_0_nuc_consensus_p8,
                  one_point_five_subs_0_nuc_consensus_p9)


percent_polymorphic = c()

#Count how many sites exist that do not have the highest mutation being a count of 100, these are polymorphic sites
for(i in site_files){
  
  data.frame <- i[,-1]
  
  highestMutation <- apply(data.frame, 1, max)
  
  polymorphic_sites = sum(highestMutation != 100 & highestMutation != 0)
  
  percent_polymorphic = c(percent_polymorphic, polymorphic_sites/length(highestMutation))
  
}

#Write out the percentage of polymorphic sites into a table
write.table(matrix(percent_polymorphic, nrow = 1), "polymorphic_percentage.csv", sep = ",", col.names = F, row.names = F)
