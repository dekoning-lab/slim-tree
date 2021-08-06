#Code adapted from Mobina Kazemi Mehrabadi
#Population (Ne), mutation rate (mu), population name (p) and output directory (outputdir) are used as arguments.
args = commandArgs(trailingOnly = TRUE)
#print(args)
Ne <- strtoi(args[1])
mu <- as.numeric(args[2])
p <- args[3]
outputdir <- args[4]


#Load codon data
codon<-read.table(paste(outputdir, "/fitnessDataFiles/slim_codon_nums.csv", sep = ""), sep = ",", header = T)
codon<-subset(codon[-c(49,51,57),],select = c(2,3))
colnames(codon)[c(1,2)]<-c("codon", "aa")

#Load Fitness profile data
Fitness_table<-as.matrix(read.delim(paste(outputdir, "/fitnessDataFiles/table_fitness_profiles.csv", sep = ""),sep = ",", header = TRUE))
fitness_initial<-matrix(data = 0,nrow = 20,ncol = 1)
rownames(fitness_initial)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

#Qmatrix of a given fitness profile
Qmatrix<-function(fitness){
  #calculating Qmatrix 2Nu*(1-exp(-sij)/1-exp(-2N*sij))
  Q_matrix<-matrix(data = 0,nrow = 61,ncol = 61)
  colnames(Q_matrix)=rownames(Q_matrix)=c("TTC","TTT","TTA","TTG","TCC","TCT","TCA","TCG","TAC","TAT","TGC","TGT","TGG","CTC","CTT","CTA","CTG","CCC","CCT" ,"CCA","CCG","CAC","CAT","CAA","CAG","CGC","CGT","CGA","CGG","ATC","ATT","ATA","ATG","ACC","ACT","ACA","ACG","AAC" ,"AAT","AAA","AAG","AGC","AGT","AGA","AGG","GTC","GTT","GTA","GTG","GCC","GCT","GCA","GCG","GAC","GAT","GAA","GAG" ,"GGC","GGT","GGA","GGG")
  for (i in rownames(Q_matrix)) {
    for (j in colnames(Q_matrix)) {
      #If the distance between two codons is 2 or 3 the Qmatrix would be zero
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==3){
        Q_matrix[i,j]<-1e-20
      }
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==2){
        Q_matrix[i,j]<-1e-20
      }
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==1){
        sij<-fitness[which(rownames(fitness_initial)==(codon$aa[which(codon$codon == j)]))]-fitness[which(rownames(fitness_initial)==(codon$aa[which(codon$codon == i)]))]
        if(sij!=0){
          Q_matrix[i,j]<-2*Ne*mu*((1-exp(-sij))/(1-exp(-2*Ne*sij)))
        }
        if(sij==0){
          Q_matrix[i,j]<-mu
        }
      }
      if(i==j){
        Q_matrix[i,j]<-0
      }
    }
  }
  #diogonal values
  for (i in rownames(Q_matrix)) {
    Q_matrix[i,i]<-(-(sum(Q_matrix[i,])))
  }
  return(Q_matrix)
}
#Stationary distribution of a given fitness profile
fitness_profile<-function(fitness){
  #Qmatrix of a given fitness profile
  ##calculating Qmatrix 2Nu*(1-exp(-sij)/1-exp(-2N*sij))
  Q_matrix<-matrix(data = 0,nrow = 61,ncol = 61)
  colnames(Q_matrix)=rownames(Q_matrix)=c("TTC","TTT","TTA","TTG","TCC","TCT","TCA","TCG","TAC","TAT","TGC","TGT","TGG","CTC","CTT","CTA","CTG","CCC","CCT" ,"CCA","CCG","CAC","CAT","CAA","CAG","CGC","CGT","CGA","CGG","ATC","ATT","ATA","ATG","ACC","ACT","ACA","ACG","AAC" ,"AAT","AAA","AAG","AGC","AGT","AGA","AGG","GTC","GTT","GTA","GTG","GCC","GCT","GCA","GCG","GAC","GAT","GAA","GAG" ,"GGC","GGT","GGA","GGG")
  for (i in rownames(Q_matrix)) {
    for (j in colnames(Q_matrix)) {
      #If the distance between two codons is 2 or 3 the Qmatrix would be zero
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==3){
        Q_matrix[i,j]<-1e-20
      }
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==2){
        Q_matrix[i,j]<-1e-20
      }
      if(length(which((unlist(strsplit(i,""))==unlist(strsplit(j,"")))==F))==1){
        sij<-fitness[which(rownames(fitness_initial)==(codon$aa[which(codon$codon == j)]))]-fitness[which(rownames(fitness_initial)==(codon$aa[which(codon$codon == i)]))]
        if(sij!=0){
          Q_matrix[i,j]<-2*Ne*mu*((1-exp(-sij))/(1-exp(-2*Ne*sij)))
        }
        if(sij==0){
          Q_matrix[i,j]<-mu
        }
      }
      if(i==j){
        Q_matrix[i,j]<-0
      }
    }
  }
  #diogonal values
  for (i in rownames(Q_matrix)) {
    Q_matrix[i,i]<-(-(sum(Q_matrix[i,])))
  }
  #Stationary distribution numerically
  coefficient_matrix<-Q_matrix
  coefficient_matrix[,61]<-1
  t_coefficient_matrix<-t(coefficient_matrix)
  Zeros<-matrix(data = 0,nrow = 61,ncol = 1)
  Zeros[61]<-1
  #Stationary distribution
  Pi<-solve(t_coefficient_matrix,Zeros,tol = 1e-40)
  Pi<-as.data.frame(Pi)
  for(i in rownames(Pi)){
    Pi[i,2]<-codon$aa[which(codon$codon==i)]
  }
  colnames(Pi)<-c("Pi","aa")
  #sum over different amino acids
  Pi_aa<-Pi %>% group_by(as.factor(aa))  %>% summarize(a=sum(Pi))
  aa<-Pi_aa$`as.factor(aa)`
  Pi_aa<-Pi_aa[,-1]
  Pi_aa[ Pi_aa < 0] = 0; # Jason
  rownames(Pi_aa)<-aa
  return(Pi_aa)
}
#amino acid stationary distribution of a given codon stationary distribution
Pi_aa<-function(Pi){
  for(i in rownames(Pi)){
    Pi[i,2]<-codon$aa[which(codon$codon==i)]
  }
  colnames(Pi)<-c("Pi","aa")
  #sum over different amino acids
  Pi_aa<-Pi %>% group_by(as.factor(aa))  %>% summarize(a=sum(Pi))
  aa<-Pi_aa$`as.factor(aa)`
  Pi_aa<-Pi_aa[,-1]
  Pi_aa[ Pi_aa < 0] = 0; # Jason
  rownames(Pi_aa)<-aa
  return(Pi_aa)
}
#codon stationary distribution of a given fitness
Stationary_distribution<-function(fitness){
  #calculating Qmatrix 2Nu*(1-exp(-sij)/1-exp(-2N*sij))
  Q_matrix<-Qmatrix(fitness)
  #Stationary distribution numerically
  coefficient_matrix<-Q_matrix
  coefficient_matrix[,61]<-1
  t_coefficient_matrix<-t(coefficient_matrix)
  Zeros<-matrix(data = 0,nrow = 61,ncol = 1)
  Zeros[61]<-1
  #Stationary distribution
  Pi<-solve(t_coefficient_matrix,Zeros,tol = 1e-40)
  Pi<-as.data.frame(Pi)
  colnames(Pi)<-"Pi"
  return(Pi)
}
#Synonymous scaling factor based on the stationary distribution and the qmatrix
scaling_factor<-function(Q_matrix,Pi){
  scaling_factor<-0
  for (i in 1:61) {
    scaling_factor<-scaling_factor+Pi[i,1]*(sum(Q_matrix[i,which(rownames(Q_matrix)%in%codon$codon[which(codon$aa==codon$aa[which(codon$codon==rownames(Q_matrix)[i])])])])-Q_matrix[i,i])
  }
  return(scaling_factor)
}
Non_synonymous_codon_finder<-function(target,Qmatrix){
  #the amino acid that our target codon belongs to
  b<-codon$aa[which(codon$codon==target)]
  #codons that are coded for the target amino acid
  c<-codon$codon[which(codon$aa==b)]
  #which codons are coded for different amino acids rather than our target one
  a<-which(!rownames(Qmatrix)%in%c)
  return(a)
}
RMSD<-function(Pi1,Pi2){
  a<-sqrt((sum((Pi1- Pi2)^2))/dim(Pi1)[1])
  return(a)
}
####Equilibrium dN/dS for different ft=itness profiles in psi-c50 file####

dN_dS_equilibrium <- matrix(0, nrow = ncol(Fitness_table), ncol = 2)

for (f in 1:(ncol(Fitness_table))){
  fitness<-Fitness_table[,f]
  #Stationary distribution and Qmatrix of the fitness profile
  Pi<-Stationary_distribution(fitness)
  Q_matrix<-Qmatrix(fitness)
  scaling_0 <- scaling_factor(Q_matrix,Pi)
  #divide the qmatrix to the scaling factor
  #Q_matrix<-Q_matrix/scaling_factor(Q_matrix,Pi)
  #initialize the a_final to zero to calculate the dN/dS equilibrium
  a_final<-0
  for (i in 1:nrow(Pi)) {
    a_final<-a_final+Pi[i,1]*sum(Q_matrix[i,Non_synonymous_codon_finder(rownames(Q_matrix)[i],Q_matrix)])
    
  }
  dN_dS_equilibrium[f, 1] <- f
  dN_dS_equilibrium[f, 2]<-a_final/scaling_0
  
}
colnames(dN_dS_equilibrium) <- c("Profile", "Value")
write.csv(dN_dS_equilibrium, paste(outputdir, "/", p , "_dNdSDistributions.csv", sep = ""), row.names = FALSE)
print(dN_dS_equilibrium)
