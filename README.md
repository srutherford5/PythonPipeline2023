# PythonPipeline2023
The first function in the code is used to create the ouput folder and also call and split all donor reads
I'm doing the differential expression track.
This pipeline requires several different imported tools as follows:
Biopython, logging, wget, sys, os, Entrez, SeqIO, pandas, and csv.
Sections of the first function calling the original SRAs have been commented out to allow
for the substitue sample file that have been provided. 
Listed just as their SRA accession numbers, these gz files should be accessible for download.
Fasta file for blast database building is also provided in the inital download and all are
moved accordingly into the output folder.
Depending on server speeds the fastq-dump is taking ~10 minutes to convert to pair-ended 
fastq files.
All code output will be written to PipelineProject_Samantha_Rutherford including the 
PipelineProject.log with all quantitative results.
