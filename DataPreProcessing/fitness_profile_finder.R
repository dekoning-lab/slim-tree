#Additional software package to find fitness profiles from stationary site-heterogeneous
#frequency profiles.

#Required arguments are -v : mutation rate, -N : population size, -f : link to site-heterogeneous
#frequency profiles csv file - amino acids should be in alphabetical order with no
#row names or col names

#Required packages:
#BB
#data.table
#dplyr


library(BB)
library(data.table)
library(dplyr)

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
args = commandArgs(trailingOnly=TRUE)
# args =  c("-v", "1e-7", "-N", "100", "-f", "Site-frequencies.csv")
cnames =  args[seq(1, length(args), 2)]
values = args[ seq(0, length(args), 2)]

args = data.frame(values)
args = transpose(args)
colnames(args) = cnames


#Get key parameters from args
v = as.numeric(args$'-v') #mutation rate
N = as.numeric(args$'-N') #population size
site_frequencies = read.csv(args$'-f', header = F)  #site heterogeneous frequency profiles


#Computes the stationary distribution of predicted amino acids and sums over different amino acid to get Pi_aa
fitness_profile<-function (fitness){
  
  #calculates Qmatrix  = 2Nv*(1-exp(-sij)/1-exp(-2N*sij))
  Q_matrix = matrix(data = 0, nrow = nCodons, ncol = nCodons)
  colnames(Q_matrix) = rownames(Q_matrix) = codons
  
  for (i in rownames(Q_matrix)) {
    
    for (j in colnames(Q_matrix)) {
      
      row_nucs = unlist(strsplit(i,""))
      col_nucs = unlist(strsplit(j,""))
      num_dif_nucs = sum(row_nucs != col_nucs)
      
      #If the two tri-nucleotide sequences differ by more then one nucleotide set Q-matrix value to 1e-20
      if(num_dif_nucs > 1){
        Q_matrix[i,j] = 1e-20 
      #If the two tri-nucleotide sequences differ by 1 difference calculate Q
      } else if (num_dif_nucs == 1){
        sij = fitness[AAs==codons_to_AAs[codons==j]] - fitness[AAs==codons_to_AAs[codons==i]]
        print(sij)
        if(sij!=0){
          Q_matrix[i,j] = 2*N*v*((1-exp(-sij))/(1-exp(-2*N*sij)))
        } else {
          Q_matrix[i,j] = v
        }
      #If the two tri-nucleotide sequences are the same - set to 0
      } else if(i==j){
        Q_matrix[i,j] = 0
      } else {
        print("This should never be reached...something is bad")
        exit()
      }
    }
  }
  
  
  #Calculate all diagonal values in the Q matrix
  for (i in rownames(Q_matrix)) {
    Q_matrix[i,i] = (-(sum(Q_matrix[i,])))
  }
  
  
  #Set up coefficient and zero matrices
  coefficient_matrix = Q_matrix
  coefficient_matrix[,nCodons] = 1
  t_coefficient_matrix = t(coefficient_matrix)
  Zeros = matrix(data = 0,nrow = nCodons,ncol = 1)
  Zeros[nCodons] = 1
  
  
  #Solve the stationary distribution
  Pi = solve(t_coefficient_matrix,Zeros,tol = 1e-40)
  Pi = as.data.frame(Pi)
  for(i in rownames(Pi)){
    Pi[i,2] = codons_to_AAs[codons==i]
  }
  colnames(Pi) = c("Pi","aa")
  
  #sum over different amino acids
  Pi_aa = Pi %>% group_by(as.factor(aa))  %>% summarize(a=sum(Pi))
  aa = Pi_aa$`as.factor(aa)`
  Pi_aa = Pi_aa[,-1]
  Pi_aa[ Pi_aa < 0] = 0;
  rownames(Pi_aa)  = AAs
  
  
  return(Pi_aa)
}

#Maximum likelihood function which finds result for a profile given its fitness
MLE <- function(fitness,profile_num){
  profile_frequencies = site_frequencies[,profile_num]
  max_frequency = which(profile_frequencies == max(profile_frequencies))
  fitness[max_frequency] = 1
  Pi_aa = fitness_profile(fitness)
  result=0
  logQ = as.data.frame(log2(Pi_aa))
  logP = as.data.frame(log2(profile_frequencies))
  colnames(logP) = "log"
  colnames(logQ) = "log"
  logP$log[!is.finite(logP$log)] = log2(1e-40)
  logQ$log[!is.finite(logQ$log)] = log2(1e-40)
  result = sum(profile_frequencies*(logP-logQ))
  return(result)
}


#Converts equilibrium frequencies to fitnesses as a starting point for optimization
equilibrium_to_fitness <- function(equilib_freq, population_size){
  pi = equilib_freq
  epsilon = 1e-316
  Ne = 2 * population_size;
  pi[equilib_freq < epsilon] = epsilon
  pi =  pi/sum(pi)
  return ((log(pi) - log(max(pi)) + Ne) / Ne);

}


#Set up final frequency and fitness tables
final_table_fitness<-c()



#Recurse through each frequency profile and find respective fitness profile
for (profile_num in 1:ncol(site_frequencies)) {
  #b is the starting point for the optimization function which is the output of equilibrium_to_fitness function
  psi = equilibrium_to_fitness2(site_frequencies[, profile_num], N)
  
  #optimization process
  if (MLE(psi,profile_num) != 0) {
    optimized_fitness = spg(psi, MLE, lower=0.1, upper=1.5, profile_num=profile_num)
    fitness_result = optimized_fitness$par
  } else {
    # Special case where starting point is exact
    fitness_result = psi;
  }
  
  # Check MSE
  diff <- (site_frequencies[,profile_num] - fitness_profile( fitness_result ))^2
  print(sum(diff)/nAAs)
  
  
  #Make combined fitness table
  final_table_fitness = cbind(final_table_fitness, fitness_result)
  print(fitness_result)
  
}


#Add stop codons at 0 frequency and lowest fitness
new_site_frequencies = rbind(site_frequencies, rep(0, ncol(site_frequencies)))
final_table_fitness = rbind(final_table_fitness, rep(min(final_table_fitness), ncol(final_table_fitness)))



#Write out the final tables
write.table(x = new_site_frequencies, file = "table_stationary_distributions.csv", col.names = F, row.names = F, sep = ",")
write.table(x = final_table_fitness, file = "table_fitness_profiles.csv", col.names = F, row.names = F, sep = ",")
