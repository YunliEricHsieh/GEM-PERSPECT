#import mean correlations
mcorr = read.csv("Results/mean_correlation_of_phynotype/mean_correlations.csv", header = TRUE)
#import mutant table
mutab = read.csv("Data/Mutant_phenotypes_table_filtered_by_num.csv", header = T)
#Remove NA corr values
mcorr=mcorr[!is.na(mcorr$Mean.Correlation),]
#find the genes with mean correlation > 0.5
geneID = mcorr$Gene[mcorr$Mean.Correlation > 0.5]
mutab = mutab[mutab$Gene %in% geneID,]

write.csv(mutab, file = "Data/Mutant_phenotypes_table_filtered_by_num_and_cor.csv", row.names = FALSE)
