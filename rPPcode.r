library(sleuth)
#Change working directory back to our output folder
setwd('PipelineProject_Samantha_Rutherford')

#Read in table describing Kallisto output
stab = read.table("Sleuth_table.txt", header=TRUE, stringsAsFactors=FALSE)
so = sleuth_prep(stab)
#Fit model to compare the two conditions
so = sleuth_fit(so, ~condition, 'full')
#Fit reduced model for likelihood ratio test
so = sleuth_fit(so, ~1, 'reduced')
#Likelihood ratio test in action
so = sleuth_lrt(so, 'reduced', 'full')

#Extract table
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
#Filter transcript at FDR<0.05
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |>dplyr::arrange(pval)
#Create text file

#Write table out to result file
write.table(sleuth_significant[,c("target_id","test_stat", "pval", "qval")], file ="fdr05_results.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
