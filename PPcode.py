#Start with imports python will need to run the code
import logging
import wget
import sys
import os
from Bio import Entrez
from Bio import SeqIO
import pandas as pd

#Part 1 Sample Accession
#Import test data by setting a stock URL variable
url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/'
#List of ids in order of patient and day 2 or 6 post infection
sra_ids = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
donors = ['Donor 1','2dpi','Donor 1','6dpi','Donor 3','2dpi','Donor 3','6dpi']

#Create function to see if output folder already exists
def out_folder(folder, ids):
    # If it already exists we do not create the folder again
    if os.path.isdir(folder)==False:
        os.mkdir(folder)
        #Change working directory to the created folder
        os.chdir('PipelineProject_Samantha_Rutherford')
        #Set up an empty list to assign urls and a for loop to create individual urls
        i=0
        sra_urls = []
        for id in ids:
            temp =  url + ids[i] + '/' + ids[i]
            sra_urls.append(temp)
            #Use wget to download sra files and os to use fastq to split the paired-end files
            wget.download(sra_urls[i])
            cmd = 'fastq-dump -I --split-files ' + ids[i]
            os.system(cmd)
            i+=1
    else:
        #Be sure output directory is the correct one
        os.chdir('PipelineProject_Samantha_Rutherford')

#Call output folder function with designated name and our ids list
out_folder('PipelineProject_Samantha_Rutherford', sra_ids)

#Create log file
logging.basicConfig(format='%(message)s', filename='PipelineProject.log', filemode='w', level='INFO')


#Step 2 Indexing

#Use Entrez to access the HCMV genome and SeqIO for later parsing
Entrez.email = 'srutherford@luc.edu'
handle = Entrez.efetch(db ='nucleotide', id = 'NC_006273.2', rettype = 'gbwithparts', retmode = 'text')
genome = next(SeqIO.parse(handle, 'genbank'))
#Create an empty list for coding sequences
cds = []
#Use for loop to access all CDS
for feature in genome.features:
    if feature.type == 'CDS':
        cds.append(feature)
#Second loop is created to remove genbank information so the sequence can be put into FASTA format
seqs = []
for feature in cds:
    seq = genome[int(feature.location.start):int(feature.location.end)]
    if feature.location.strand == -1:
        seq.seq = seq.seq.reverse_complement()
    seq.description = feature.qualifiers['product'][0]
    seqs.append(seq)
#Write coding sequences to a fasta file
SeqIO.write(seqs, 'NC_006273.2.fasta', 'fasta')
#Create variable
logging.info('The HCMV genome (NC_006273.2) has %s CDS.', len(seqs))

#Function to ensure index and abundance files are only created once
def quant_out(folder, ids):
    #If it already exists we do not create the folder again
    if os.path.isdir(folder)==False:
        os.mkdir(folder)
        #Build the HCMV index for Kallisto with coding sequence fasta
        os.system('time kallisto index -i NC_006273.2.idx NC_006273.2.fasta --make-unique -k 31')
        #i is my favorite counter
        i = 0
        #For loop iterates through each sample, creating TPM output
        for id in ids:
            fastq1 = ids[i] + '_1.fastq'
            fastq2 = ids[i] + '_2.fastq'
            cmd = 'time kallisto quant -i NC_006273.2.idx -o results/' + ids[i] + ' -b 30 -t 2 ' + fastq1 + ' ' + fastq2
            os.system(cmd)
            i+=1

quant_out('results', sra_ids)

#Function to calculate mean, median, min, and max for each sample and write results to log
def quant_results(ids, sample):
    #Write log header and column values
    logging.info("Quantitative Statistics of HCMV Sample's TPM")
    mess = 'sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm'
    logging.info(mess)
    i = 0
    d = 0
    #Assign file variable outside of loop because each file has the same name
    call_file = 'abundance.tsv'
    #Change currnet directory the created results folder
    os.chdir('results')
    #For loop opens sample tsv results and calculates requested data
    for id in ids:
        #Change directory to respective sample and read in file with pandas
        os.chdir(ids[i])
        file = pd.read_table(call_file)
        #Specify tpm column and calculate stats
        tpm = file['tpm']
        tmin = tpm.min()
        tmed = tpm.median()
        tmean = tpm.mean()
        tmax = tpm.max()
        #Log stats and move directory back to results before next iteration
        log_mess = sample[d] + '\t' + sample[d+1] + '\t' + str(tmin) +  '\t' + str(tmed) + '\t' + str(tmean) + '\t' + str(tmax)
        logging.info(log_mess)
        os.chdir('..')
        i+=1
        d+=2

quant_results(sra_ids, donors)

#Change directory to output to create kallisto output table
os.chdir('..')
#Pandas dataframe created to ensure delimited file and written out to txt file
kdf = pd.DataFrame({'sample': sra_ids, 'condition':['2dpi','6dpi','2dpi','6dpi'], 'path': ['results/SRR5660030', 'results/SRR5660033', 'results/SRR5660044', 'results/SRR5660045']})
kout = 'Sleuth_table.txt'
with open(kout, 'w') as f:
    kstring = kdf.to_string(header=True, index=False)
    f.write(kstring)
#Change directory to call the Rscript for sleuth
os.chdir('..')
os.system('Rscript rPPcode.r')
#Switch back to output so we can log the output from sleuth
os.chdir('PipelineProject_Samantha_Rutherford')
open_sleuth = open('fdr05_results.txt', 'r')
sleuth = open_sleuth.read()
logging.info('\n')
logging.info('Differential Gene Expression Between Two Timepoints 2pi and 6pi')
logging.info(sleuth)

