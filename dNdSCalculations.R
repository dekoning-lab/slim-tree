#library(data.table)

#Set up amino acids and translation from codons to amino acids
AAs = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
codons = c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT", "TAC", "TGT", 
           "TGC", "TGG","CTT", "CTC", "CTA","CTG", "CCT", "CCC","CCA", "CCG", "CAT", "CAC",
           "CAA","CAG","CGT", "CGC","CGA","CGG","ATT", "ATC","ATA","ATG","ACT","ACC","ACA",
           "ACG","AAT","AAC" ,"AAA","AAG", "AGT", "AGC", "AGA","AGG", "GTT", "GTC","GTA","GTG",
           "GCT", "GCC","GCA","GCG", "GAT", "GAC","GAA", "GAG", "GGT", "GGC","GGA","GGG")
codons_to_AAs = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "C", "C",
                  "W","L", "L", "L", "L", "P", "P", "P", "P","H", "H", "Q", "Q", "R", "R", "R", 
                  "R","I", "I", "I", "M", "T", "T", "T",  "T", "N", "N", "K", "K", "S", "S", 
                  "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", 
                  "G", "G", "G", "G")


#Set up initial fitness matrix
nAAs = 20
nCodons = 61
fitness_mat<-matrix(data = 0,nrow = nAAs,ncol = 1)
rownames(fitness_mat)<-AAs



#Read in user input of mutation rate, population size and site-heterogeneous frequency profiles
args = commandArgs(trailingOnly = TRUE)
N <- strtoi(args[1]) #Population size
v <- as.numeric(args[2]) #Mutation rate
p <- args[3] #Population name (for output)
inputdir <- args[4] #Input directory
outputdir <- args[5] #Output directory

#Read in frequency and fitness profiles from designated locations, ignore stop codons
frequency_profiles = read.csv(paste(inputdir, "/fitnessDataFiles/table_stationary_distributions.csv", sep = ""), header = F)  #site heterogeneous frequency profiles
frequency_profiles <- frequency_profiles[1:20,]
fitness_profiles = read.csv(paste(inputdir, "/fitnessDataFiles/table_fitness_profiles.csv", sep = ""), header = F)  #site heterogeneous fitness profiles
fitness_profiles <- fitness_profiles[1:20,]

find_matrices <- function(fitness){
  #calculates Qmatrix  = 2Nv*(1-exp(-sij)/1-exp(-2N*sij))
  
  Q_matrix = matrix(data = 0, nrow = nCodons, ncol = nCodons)
  colnames(Q_matrix) = rownames(Q_matrix) = codons
  
  V_matrix = matrix(data = 0, nrow = nCodons, ncol = nCodons)
  colnames(V_matrix) = rownames(V_matrix) = codons
  
  for (i in rownames(Q_matrix)) {
    
    for (j in colnames(Q_matrix)) {
      
      row_nucs = unlist(strsplit(i,""))
      col_nucs = unlist(strsplit(j,""))
      num_dif_nucs = sum(row_nucs != col_nucs)
      
      V_matrix[i,j] = (v)
      
      
      #If the two tri-nucleotide sequences differ by more then one nucleotide set Q-matrix value to 1e-20
      if(num_dif_nucs > 1){
        Q_matrix[i,j] = 1e-20 
        #If the two tri-nucleotide sequences differ by 1 difference calculate Q
      } else if (num_dif_nucs == 1){
        sij = fitness[AAs==codons_to_AAs[codons==j]] - fitness[AAs==codons_to_AAs[codons==i]]
        if(sij!=0){
          Q_matrix[i,j] = 2*N*v*((1-exp(-sij))/(1-exp(-2*N*sij)))
        } else {
          Q_matrix[i,j] = v
        }
        #If the two tri-nucleotide sequences are the same - set to 0
      } else if(i==j){
        Q_matrix[i,j] = 0
        V_matrix[i,j] = 0
      } else {
        print("This should never be reached...something is bad")
        exit()
      }
    }
  }
  
  
  
  return_list = list("Q_matrix" = Q_matrix, "v_matrix" = V_matrix)
  return(return_list)
}


dnds_vals = c()

#Recurse through each fitness profile and find dnds value
for (profile_num in 1:ncol(fitness_profiles)) {
  
  #Get current profile's fitnesses and frequencies
  current_fitnesses = unlist(fitness_profiles[profile_num])
  names(current_fitnesses) = AAs
  current_frequencies = unlist(frequency_profiles[profile_num])
  names(current_frequencies) = AAs
  
  #Get Q and V matrices
  matrices = find_matrices(current_fitnesses)
  pn = c()
  ps = c()
  
  p_star_n = c()
  p_star_s = c()
  
  #Recurse through each AA and find pn, ps, pn* and ps*
  for (aa in AAs){
    available_codons = (codons_to_AAs == aa) #Determines all synonymous and non-synonymous codons for that amino acid
    
    n_codons = sum(available_codons)
    frequency = current_frequencies[aa] #Get frequency of that amino acid in the profile
    
    if(n_codons > 1){
      ps = c(ps, frequency*rowSums(matrices$Q_matrix[available_codons,available_codons]))
      p_star_s = c(p_star_s, frequency*rowSums(matrices$v_matrix[available_codons,available_codons]))
    }
    
    pn = c(pn, frequency*rowSums(matrices$Q_matrix[!available_codons,!available_codons]))
    p_star_n = c(p_star_n, frequency*rowSums(matrices$v_matrix[!available_codons,!available_codons]))
    
    
  }
  
  #Calculate dN/dS
  dnds = (sum(pn)*sum(p_star_s))/(sum(ps)*sum(p_star_n))
  dnds_vals = c(dnds_vals, dnds)
}

#Write output in a csv file
write.csv(dnds_vals, paste(outputdir, "/", p , "_dNdSDistributions.csv", sep = ""))
