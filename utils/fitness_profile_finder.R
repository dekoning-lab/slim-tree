#Finds fitness profiles from stationary site-heterogeneous frequency profiles.

#Required arguments are -v : mutation rate, -N : population size, -f : link to site-heterogeneous
#frequency profiles csv file - amino acids should be in alphabetical order with no
#row names or col names

#Required packages:
#dplyr
#BB
#data.table
#optparse
#seqinr


#Install packages not already installed
list.of.packages <- c("dplyr", "BB", "data.table", "optparse", "seqinr", "doParallel", "Rfast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

#load packages
suppressMessages(library(dplyr))
suppressMessages(library(BB))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(seqinr))
suppressMessages(library(doParallel))
suppressMessages(library(Rfast))


#Read in user input of mutation rate, population size and site-heterogeneous frequency profiles
list_of_options = list(
  make_option(c("-N", "--population_size"), type = "integer",
              help = "Population size of the simulation"),
  make_option(c("-v", "--mutation_rate"), type = "numeric",
              help = "Mutation rate for the simulation"),
  make_option(c("-f", "--codon_frequencies_file"), type = "character",
              help = "File containing site-heterogeneous frequency profiles"),
  make_option(c("-o", "--output_file"), type = "character",
              help = "Output file for final fitness profiles"))

input_parser = OptionParser(option_list = list_of_options)
args = parse_args(input_parser)
v = args$mutation_rate
N = args$population_size


#Function to find the frequency of each amino acid
find_amino_acid_freq <- function(AA, codon_frequencies){
  freqs = codon_frequencies[codon_frequencies$AA == AA,]
  freqs = colSums(freqs[,-ncol(freqs)])
  return(freqs)
}


#Function to convert codons to an amino acid
codon_to_AA <- function(codon_name){
  #turn codon name to a list of characters in codon
  codon_name <- tolower(codon_name)
  codon_list <- strsplit(codon_name, split = "")[[1]]

  #if codons given with u instead of t, convert to t
  codon_list[codon_list == "u"] <- "t"

  #translate to amino acids
  AA = translate(codon_list)
  
  return(AA)
}

#Function to convert list of codons to amino acids
codons_to_AAs <- function(codons){
  return(sapply(codons, codon_to_AA))
}


#Open site frequencies file and convert to amino acid frequencies
codon_frequencies = read.csv(args$codon_frequencies_file, header = F, row.names = 1)
codon_frequencies$AA = codons_to_AAs(rownames(codon_frequencies))

#Get AAs and codons in frequency file
AAs = unique(codon_frequencies$AA)
codons = rownames(codon_frequencies)

nAAs = length(AAs)
nCodons = length(codons)


#Find amino acid frequencies
aa_frequencies = t(sapply(AAs, find_amino_acid_freq, codon_frequencies = codon_frequencies))



#Set up initial fitness matrix
fitness_mat<-matrix(data = 0,nrow = nAAs,ncol = 1)
rownames(fitness_mat)<-AAs



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
        sij = fitness[AAs==codons_to_AAs(j)] - fitness[AAs==codons_to_AAs(i)]
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
    Pi[i,2] = codons_to_AAs(i)
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
  profile_frequencies = aa_frequencies[,profile_num]
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


registerDoParallel(cores=detectCores() - 3)

#Recurse through each frequency profile in parallel and find respective fitness profile
final_table_fitness<- foreach(profile_num =  1:ncol(aa_frequencies), 
                              .combine=cbind, .packages = c('seqinr', 'dplyr', 'BB'))  %dopar% {
                                
  #b is the starting point for the optimization function which is the output of equilibrium_to_fitness function
  psi = equilibrium_to_fitness(aa_frequencies[, profile_num], N)
  
  #optimization process
  if (MLE(psi,profile_num) != 0) {
    optimized_fitness = spg(psi, MLE, lower=0.1, upper=1.5, profile_num=profile_num)
    fitness_result = optimized_fitness$par
  } else {
    # Special case where starting point is exact
    fitness_result = psi;
  }

  #Make combined fitness table
  return(fitness_result)

}


#Set stop codon to the lowest fitness value
final_table_fitness <- rbind(final_table_fitness, colMins(final_table_fitness, value = T))
rownames(final_table_fitness) <- c(AAs, "X")



#Write out the fitness table
write.table(x = final_table_fitness, file = args$output_file, col.names = F, sep = ",")







