library(pracma)
library(mosaic)
library(resample)
library(gridExtra)
library(readxl)


#Create rows to format data for an output table
blank_line <- c("","","","","","","", "")
title_0.001 <- c("?? = 0.001","","","","","","", "")
title_0.05 <- c("?? = 0.05","","","","","","", "")
title_0.1 <- c("?? = 0.1","","","","","","", "")
title_0.2 <- c("?? = 0.2","","","","","","", "")



#Obtain the branch lengths for each dataset - note these branch lengths were summarized in excel
summarized_branches_0.001 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                        sheet = "0.001", col_types = c("skip", "numeric", "numeric", 
                                                                       "numeric", "numeric", "numeric", 
                                                                       "numeric", "numeric", "numeric"))
summarized_branches_0.001$p9 <- summarized_branches_0.001$p9 + summarized_branches_0.001$p2
summarized_branches_0.001 <- summarized_branches_0.001[ , order(colnames(summarized_branches_0.001))]
summarized_branches_0.05 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                       sheet = "0.05", col_types = c("skip", "numeric", "numeric", 
                                                                     "numeric", "numeric", "numeric", 
                                                                     "numeric", "numeric", "numeric"))
summarized_branches_0.05$p9 <- summarized_branches_0.05$p9 + summarized_branches_0.05$p2
summarized_branches_0.05 <-summarized_branches_0.05[ , order(colnames(summarized_branches_0.05))]
summarized_branches_0.1 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                      sheet = "0.1", col_types = c("skip", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric"))
summarized_branches_0.1$p9 <- summarized_branches_0.1$p9 + summarized_branches_0.1$p2
summarized_branches_0.1 <-summarized_branches_0.1[ , order(colnames(summarized_branches_0.1))]
summarized_branches_0.2 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                      sheet = "0.2", col_types = c("skip", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric"))
summarized_branches_0.2$p9 <- summarized_branches_0.2$p9 + summarized_branches_0.2$p2
summarized_branches_0.2 <-summarized_branches_0.2[ , order(colnames(summarized_branches_0.2))]
summarized_branches_kimura_0.001 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                               sheet = "k0.001", col_types = c("skip", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric", 
                                                                               "numeric", "numeric", "numeric"))
summarized_branches_kimura_0.001$p9 <- summarized_branches_kimura_0.001$p9 + summarized_branches_kimura_0.001$p2
summarized_branches_kimura_0.001 <-summarized_branches_kimura_0.001[ , order(colnames(summarized_branches_kimura_0.001))]
summarized_branches_kimura_0.05 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                              sheet = "k0.05", col_types = c("skip", "numeric", "numeric", 
                                                                             "numeric", "numeric", "numeric", 
                                                                             "numeric", "numeric", "numeric"))
summarized_branches_kimura_0.05$p9 <- summarized_branches_kimura_0.05$p9 + summarized_branches_kimura_0.05$p2
summarized_branches_kimura_0.05 <-summarized_branches_kimura_0.05[ , order(colnames(summarized_branches_kimura_0.05))]
summarized_branches_kimura_0.1 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                             sheet = "k0.1", col_types = c("skip", "numeric", "numeric", 
                                                                           "numeric", "numeric", "numeric", 
                                                                           "numeric", "numeric", "numeric"))
summarized_branches_kimura_0.1$p9 <- summarized_branches_kimura_0.1$p9 + summarized_branches_kimura_0.1$p2
summarized_branches_kimura_0.1 <-summarized_branches_kimura_0.1[ , order(colnames(summarized_branches_kimura_0.1))]
summarized_branches_kimura_0.2 <- read_excel("C:/Users/ebrin/Dropbox (APJdKL)/MDSC_508/PostProcessing/branchlengthsGood/summarized_branches.xlsx", 
                                             sheet = "k0.2", col_types = c("skip", "numeric", "numeric", 
                                                                           "numeric", "numeric", "numeric", 
                                                                           "numeric", "numeric", "numeric"))

summarized_branches_kimura_0.2$p9 <- summarized_branches_kimura_0.2$p9 + summarized_branches_kimura_0.2$p2
summarized_branches_kimura_0.2 <-summarized_branches_kimura_0.2[ , order(colnames(summarized_branches_kimura_0.2))]                                                                   

#Run tests for normality of data

runNormalityTests <- function(pop){
  pop_3 <- shapiro.test(pop$p3)
  pop_4 <- shapiro.test(pop$p4)
  pop_5 <- shapiro.test(pop$p5)
  pop_6 <- shapiro.test(pop$p6)
  pop_7 <- shapiro.test(pop$p7)
  pop_8 <- shapiro.test(pop$p8)
  pop_9 <- shapiro.test(pop$p9)
  
  p_values <- c(pop_3$p.value,pop_4$p.value, pop_5$p.value, pop_6$p.value, pop_7$p.value, pop_8$p.value, pop_9$p.value)
  statistics <- c(pop_3$statistic,pop_4$statistic, pop_5$statistic, pop_6$statistic, pop_7$statistic, pop_8$statistic, pop_9$statistic)
  
  tests <- rbind(p_values,statistics)
  
  return(tests)
}

normality_tests_kimura_0.001 <-runNormalityTests(summarized_branches_kimura_0.001)
normality_tests_kimura_0.05 <- runNormalityTests(summarized_branches_kimura_0.05)
normality_tests_kimura_0.1 <-  runNormalityTests(summarized_branches_kimura_0.1)
normality_tests_kimura_0.2 <-  runNormalityTests(summarized_branches_kimura_0.2)

formatted_normality_kimura <- rbind(blank_line,title_0.001, normality_tests_kimura_0.001, blank_line, 
                                    title_0.05, normality_tests_kimura_0.05, blank_line,
                                    title_0.1, normality_tests_kimura_0.1, blank_line,
                                    title_0.2, normality_tests_kimura_0.2)
colnames(formatted_normality_kimura) <- c("p3", "p4", "p5", "p6", "p7", "p8", "p2/9")

write.csv( formatted_normality_kimura, file = "normality_test_branches_kimura.csv")


normality_tests_0.001 <-runNormalityTests(summarized_branches_0.001)
normality_tests_0.05 <- runNormalityTests(summarized_branches_0.05)
normality_tests_0.1 <-  runNormalityTests(summarized_branches_0.1)
normality_tests_0.2 <-  runNormalityTests(summarized_branches_0.2)

formatted_normality<- rbind(blank_line,title_0.001, normality_tests_0.001, blank_line, 
                            title_0.05, normality_tests_0.05, blank_line,
                            title_0.1, normality_tests_0.1, blank_line,
                            title_0.2, normality_tests_0.2)
colnames(formatted_normality) <- c("p3", "p4", "p5", "p6", "p7", "p8", "p2/9")

write.csv( formatted_normality, file = "normality_test_branches.csv")



#Expected Branch Lengths                                                                                                                                                                                                    #Figure out the expected length of each branch
p_3 <- 0.094
p_4 <- 0.0384
p_5 <- 0.0576 
p_6 <- 0.0576 
p_7 <- 0.096 
p_8 <- 0.19 
p_9 <- 0.41 
total <- sum(c(p_3, p_4, p_5, p_6, p_7, p_8, p_9))

expected_means <- c(p_3, p_4, p_5, p_6, p_7, p_8, p_9, total)


#Make Table of Statistical Significance of lengths
runTests <- function(pop){
  pop <- pop[,-1]
  pop$totals <- rowSums(pop) 
  pop_3 <- wilcox.test(pop$p3, mu = p_3, alternative = "two.sided")
  pop_4 <- wilcox.test(pop$p4, mu = p_4, alternative = "two.sided")
  pop_5 <- wilcox.test(pop$p5, mu = p_5, alternative = "two.sided")
  pop_6 <- wilcox.test(pop$p6, mu = p_6, alternative = "two.sided")
  pop_7 <- wilcox.test(pop$p7, mu = p_7, alternative = "two.sided")
  pop_8 <- wilcox.test(pop$p8, mu = p_8, alternative = "two.sided")
  pop_9 <- wilcox.test(pop$p9, mu = p_9, alternative = "two.sided")
  all_pop <- wilcox.test(pop$totals, mu = total, alternative = "two.sided")
  
  tests <- cbind(pop_3,pop_4, pop_5, pop_6, pop_7, pop_8, pop_9, all_pop)
  
  return(tests)
}


#Make a table of statistics that include the Wilcox test and error values
extractRequiredStatistics <- function(pop_dat, w_test_dat){
  pop_dat <- pop_dat[,-1]
  pop_dat$totals <- rowSums(pop_dat)
  
  percent_error_values_col <- ((pop_dat[,1] - expected_means[1])/expected_means[1]) * 100
  percent_error_values <- paste(signif(colMeans(percent_error_values_col), digits = 3), " � ",
                                signif(colStdevs(percent_error_values_col), digits = 3))
  
  for (i in 2:ncol(pop_dat)){
    percent_error_values_col <- ((pop_dat[,i] - expected_means[i])/expected_means[i]) * 100
    percent_error_values <- cbind(percent_error_values, paste(signif(colMeans(percent_error_values_col), digits = 3), "�",
                                                              signif(colStdevs(percent_error_values_col), digits = 3)))
  }
  
  
  means <- paste(signif(colMeans(pop_dat),digits=4))
  sds <- paste(signif(colStdevs(pop_dat),digits=4))
  w_stats <- signif(unlist(w_test_dat[1,]), digits = 3)
  p_values <- signif(unlist(w_test_dat[3,]), digits = 3)
  
  formatted_data <- rbind(means, sds, w_stats, p_values, percent_error_values)#, full_chi)
  rownames(formatted_data) <- c("Mean (n = 100)", "Standard Deviation", "Statistic (W)", "p value", "MAPE(%)")#, "")
  colnames(formatted_data) <- c("p3", "p4", "p5", "p6", "p7", "p8", "p2/9", "All Branches")
  
  return(formatted_data)
  
}





#Get statistics for the kimura data and make into a table
tests_kimura_0.001 <- extractRequiredStatistics(summarized_branches_kimura_0.001 , runTests(summarized_branches_kimura_0.001))
tests_kimura_0.05 <- extractRequiredStatistics(summarized_branches_kimura_0.05, runTests(summarized_branches_kimura_0.05))
tests_kimura_0.1 <- extractRequiredStatistics(summarized_branches_kimura_0.1, runTests(summarized_branches_kimura_0.1))
tests_kimura_0.2 <- extractRequiredStatistics(summarized_branches_kimura_0.2, runTests(summarized_branches_kimura_0.2))

test_statistics_kimura <- rbind(expected_means, blank_line, title_0.001, tests_kimura_0.001,
                                  blank_line, title_0.05, tests_kimura_0.05, blank_line,
                                  title_0.1, tests_kimura_0.1, blank_line, title_0.2, tests_kimura_0.2)
write.csv(test_statistics_kimura, file = "test statistics kimura.csv")


#Get statistics for the non-kimura data and make into a table
tests_0.001 <- extractRequiredStatistics(summarized_branches_0.001 , runTests(summarized_branches_0.001))
tests_0.05 <- extractRequiredStatistics(summarized_branches_0.05, runTests(summarized_branches_0.05))
tests_0.1 <- extractRequiredStatistics(summarized_branches_0.1, runTests(summarized_branches_0.1))
tests_0.2 <- extractRequiredStatistics(summarized_branches_0.2, runTests(summarized_branches_0.2))

test_statistics <- rbind(expected_means, blank_line, title_0.001, tests_0.001,
                           blank_line, title_0.05, tests_0.05, blank_line,
                           title_0.1, tests_0.1, blank_line, title_0.2, tests_0.2)
write.csv(test_statistics, file = "test statistics.csv")





#Find the error values 
find_error <- function(pop_dat){
  pop_dat <- pop_dat[,-1]
  percent_error_values_col <- ((pop_dat[,1] - expected_means[1])/expected_means[1]) * 100
  percent_error_values <- c(unlist(percent_error_values_col))
  
  for (i in 2:ncol(pop_dat)){
    percent_error_values_col <- ((pop_dat[,i] - expected_means[i])/expected_means[i]) * 100
    percent_error_values <- c(percent_error_values, unlist(percent_error_values_col))
  }
  return(percent_error_values)
}


substitutions <- c(rep(p_3, 100), rep(p_4, 100), rep(p_5, 100), rep(p_6, 100), rep(p_7, 100), rep(p_8, 100), rep(p_9, 100))



#Get the kimura error and format into data frame
kimura_0.001_error <- find_error(summarized_branches_kimura_0.001)
kimura_0.05_error <- find_error(summarized_branches_kimura_0.05)
kimura_0.1_error <- find_error(summarized_branches_kimura_0.1)
kimura_0.2_error <- find_error(summarized_branches_kimura_0.2)

kimura_percent_error_values <- data.frame(cbind(substitutions, kimura_0.001_error, kimura_0.05_error, kimura_0.1_error,
                                                kimura_0.2_error))

#Get the non-kimura error and format into data frame
not_kimura_0.001_error <- find_error(summarized_branches_0.001)
not_kimura_0.05_error <- find_error(summarized_branches_0.05)
not_kimura_0.1_error <- find_error(summarized_branches_0.1)
not_kimura_0.2_error <- find_error(summarized_branches_0.2)


not_kimura_percent_error_values <- data.frame(cbind(substitutions, not_kimura_0.001_error, not_kimura_0.05_error, 
                                                    not_kimura_0.1_error, not_kimura_0.2_error))


#Plot data in terms of branch length 
colnames(kimura_percent_error_values) <- c("Expected_Substitutions", "Error_0.001","Error_0.05" ,"Error_0.1", "Error_0.2")


plot1 <- ggplot(data = kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.001)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("A") + ylim(-110, 2100)
plot2 <- ggplot(data = kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.05)) + geom_boxplot() +xlab("Expected Branche Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("B")+ ylim(-110, 2100)
plot3 <- ggplot(data = kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.1)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("C")+ ylim(-110, 2100)
plot4 <- ggplot(data = kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.2)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("D")+ ylim(-110,2100)

grid.arrange(plot1,plot2,plot3, plot4, nrow = 2)

colnames(not_kimura_percent_error_values) <- c("Expected_Substitutions", "Error_0.001","Error_0.05" ,"Error_0.1", "Error_0.2")

plot1 <- ggplot(data = not_kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.001)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("A") + ylim(-110, 2050)
plot2 <- ggplot(data = not_kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.05)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)") + ylab("Percent Error (%)")+ggtitle("B")+ ylim(-110, 2050)
plot3 <- ggplot(data = not_kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.1)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)")+ ylab("Percent Error (%)")+ggtitle("C")+ ylim(-110, 2050)
plot4 <- ggplot(data = not_kimura_percent_error_values, aes(x=factor(Expected_Substitutions), y =  Error_0.2)) + geom_boxplot() +xlab("Expected Branch Length (subs/site)")+ ylab("Percent Error (%)")+ggtitle("D")+ ylim(-110,2050)

grid.arrange(plot1,plot2,plot3, plot4, nrow = 2)



#Make plots comparing error values of data
kimura_percent_error_values_combined <- data.frame(unlist(kimura_percent_error_values[,-1]))
colnames(kimura_percent_error_values_combined) <- c("percent")
kimura_percent_error_values_combined$absolute <- abs(kimura_percent_error_values_combined$percent)
kimura_percent_error_values_combined$names <- c(rep("0.001", 700), rep("0.05", 700), rep("0.1", 700), rep("0.2", 700))

percent_error_values_combined <- data.frame(unlist(not_kimura_percent_error_values[,-1]))
colnames(percent_error_values_combined) <- c("percent")
percent_error_values_combined$absolute <- abs(percent_error_values_combined$percent)
percent_error_values_combined$names <- c(rep("0.001", 700), rep("0.05", 700), rep("0.1", 700), rep("0.2", 700))


plot_1_kimura <- ggplot(data = kimura_percent_error_values_combined, aes(x = names, y =  percent)) + geom_boxplot() +ylab("Percent Error (%)") + xlab("Population Mutation Rate Parameter (??)")+ggtitle("A") + ylim(-110, 2100)
plot_2_kimura <- ggplot(data = kimura_percent_error_values_combined, aes(x = names, y =  absolute)) + geom_boxplot() +ylab("Absolute Percent Error (%)") + xlab("Population Mutation Rate Parameter (??)")+ggtitle("B") + ylim(-110, 2100)

grid.arrange(plot_1_kimura, plot_2_kimura, nrow = 1)

plot_1_not_kimura <- ggplot(data = percent_error_values_combined, aes(x = names, y =  percent)) + geom_boxplot() +ylab("Percent Error (%)") + xlab("Population Mutation Rate Parameter (??)")+ggtitle("A") 
plot_2_not_kimura <- ggplot(data = percent_error_values_combined, aes(x = names, y =  absolute)) + geom_boxplot() +ylab("Absolute Percent Error (%)") + xlab("Population Mutation Rate Parameter (??)")+ggtitle("B") 

grid.arrange(plot_1_not_kimura, plot_2_not_kimura, nrow = 1)


#Compare error values in a statistically significant manner
wilcox.test(abs(not_kimura_0.001_error), abs(kimura_0.001_error), alternative = "greater")
wilcox.test(abs(not_kimura_0.05_error), abs(kimura_0.05_error), alternative = "greater")
wilcox.test(abs(not_kimura_0.1_error), abs(kimura_0.1_error), alternative = "greater")
wilcox.test(abs(not_kimura_0.2_error), abs(kimura_0.2_error), alternative = "greater")


wilcox.test(not_kimura_0.001_error, kimura_0.001_error, alternative = "greater")
wilcox.test(not_kimura_0.05_error, kimura_0.05_error, alternative = "greater")
wilcox.test(not_kimura_0.1_error, kimura_0.1_error, alternative = "greater")
wilcox.test(not_kimura_0.2_error, kimura_0.2_error, alternative = "greater")


#Compare error values between theta

wilcox.test(not_kimura_0.001_error, not_kimura_0.05_error, alternative = "less")
wilcox.test(not_kimura_0.05_error, not_kimura_0.1_error, alternative = "less")
wilcox.test(not_kimura_0.1_error, not_kimura_0.2_error, alternative = "less")


wilcox.test(kimura_0.001_error, kimura_0.05_error, alternative = "less")
wilcox.test(kimura_0.05_error, kimura_0.1_error, alternative = "less")
wilcox.test(kimura_0.1_error, kimura_0.2_error, alternative = "less")

